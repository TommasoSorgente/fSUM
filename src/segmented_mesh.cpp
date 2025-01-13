#include "segmented_mesh.h"

/**********************************************************************/
/*************************** PUBLIC METHODS ***************************/
/**********************************************************************/

template<class M, class V,  class E, class P> inline
void SegmentedMesh::init(AbstractMesh<M,V,E,P> &m, Parameters &Par)
{
    prof.push("FESA::init");
    verbose = Par.get_VERBOSE();

    // load mesh and scalar fields
    m.load(Par.get_MESH_PATH().c_str());

    load_field(m, Par.get_FIELD_PATH(), Par.get_FIELD_TYPE());
    if (Par.get_FIELD_TYPE() == 1)
        extend_field_to_verts(m);
    else if (Par.get_FIELD_TYPE() == 2)
        extend_field_to_cells(m);
    else
        assert(false && "unknown field type");

    if (Par.get_FGLOBAL()) {
        load_global_field(Par.get_FIELDG_PATH());
    } else {
        for (uint i=0; i<field.size(); ++i) {
            global_field.push_back(field[i]);
        }
        auto minmax = minmax_element(global_field.begin(), global_field.end());
        field_min = *minmax.first;
        field_max = *minmax.second;
    }
    if (Par.get_GUI()) {
        color_mesh(m);
    }
    // load_second_field(m, Par.get_FIELD2_PATH());

    // create and open the output directory
    size_t last = Par.get_MESH_PATH().find_last_of('/');
    std::string name_tmp = Par.get_MESH_PATH().substr(0, last);
    size_t first = name_tmp.find_last_of('/');
    std::string name = name_tmp.substr(first+1);
    output_path = std::string(HOME_PATH) + Par.get_OUT_PATH() + name;
    open_directory(output_path);

    // print the initial field
    std::ofstream fp(output_path + "/input_field.csv");
    assert(fp.is_open());
    fp << std::fixed;
    fp << "# x, y, z, field" << std::endl;
    for (uint pid=0; pid<m.num_polys(); ++pid) {
        vec3d c = m.poly_centroid(pid);
        fp.precision(2);
        fp << c.x() << ", " << c.y() << ", " << c.z() << ", ";
        fp.precision(6);
        fp << field[pid] << std::endl;
    }
    fp.close();

    prof.pop();
}

/**********************************************************************/

template<class M, class V,  class E, class P> inline
void SegmentedMesh::segment(AbstractMesh<M,V,E,P> &m, Parameters &Par)
{
    prof.push("FESA::segment");
    compute_isovalues(Par.get_N_REGIONS(), Par.get_ISOVAL_TYPE(), Par.get_ISOVAL_VALS());

    compute_isoregions(m, Par.get_DENOISE());

    separate_ccs(m); // from now on, we switch from regions to subregions!

    prof.pop();
}

/**********************************************************************/

template<class M, class V,  class E, class P> inline
void SegmentedMesh::clean(AbstractMesh<M,V,E,P> &m, Parameters &Par)
{
    prof.push("FESA::clean");
    int i = 0, n_labels_tmp = 0;
    int n_labels_init = polys_in_subregion.size();

    // compute mesh total mass (area or volume)
    double mass = 0;
    for (uint pid=0; pid<m.num_polys(); ++pid) {
        mass += m.poly_mass(pid);
    }
    double CLEAN_THRESH = mass * Par.get_CLEAN_THRESH() / 100.;

    while (n_labels_tmp != (int)polys_in_subregion.size() && i < Par.get_N_ITER()) {
        n_labels_tmp = polys_in_subregion.size();
        remove_small_labels(m, CLEAN_THRESH);
        ++i;
    }

    if (n_labels_tmp == n_labels_init) {
        cout << "FESA Warning: did not remove any label. Filter_thresh is too high?";
    }
    prof.pop(true, "\tfound " + std::to_string(polys_in_subregion.size()) + " components");
}

/**********************************************************************/

template<class M, class V,  class E, class P> inline
void SegmentedMesh::smooth(AbstractMesh<M,V,E,P> &m, Parameters &Par)
{
    prof.push("FESA::smooth boundaries");
    int j = 0, count = 1;
    while (count > 0 && j < Par.get_N_ITER()) {
        count = smooth_boundaries(m);
        ++j;
    }
    prof.pop();
}

/**********************************************************************/

template<class M, class V,  class E, class P> inline
void SegmentedMesh::output(Polygonmesh<M,V,E,P> &m, Parameters &Par)
{
    prof.push("FESA::output");
    update_regions_map(m, polys_in_subregion);

    // print field and isovalues info
    print_field_csv(m, output_path + "/output_field.csv", prof);
    print_isovals_csv(m, percentiles, isovals, output_path + "/isovalues.txt", prof);

    // save mesh
    std::string domain_dir = output_path + "/domain";
    open_directory(domain_dir);
    m.save((domain_dir + "/domain.obj").c_str());

    // print verts and cells data (label and field) in a csv file
    print_verts_csv(m, domain_dir + "/domain_verts_data.csv", prof);
    print_cells_csv(m, subregion_region_map, domain_dir + "/domain_cells_data.csv", prof);

    // print general statistics
    Statistics stats;
    stats.init(percentiles, isovals, prof);
    stats.compute_stats(m, polys_in_subregion);

    // generate a shapefile containing the iso-regions boundaries
    print_regions_shp(m, stats, polys_in_subregion, domain_dir + "/domain", prof);

    if (Par.get_OUT_LEVEL() > 1) {
        if (Par.get_OUT_LEVEL() > 2) {
            // save a mesh for each iso-subregion
            stats.init(percentiles, isovals, prof);
            print_regions(m, subregion_region_map, polys_in_subregion, stats, output_path + "/subregions", prof);
        }
        // save a mesh for each iso-region
        restore_original_labels(m);
        print_regions(m, subregion_region_map, polys_in_region, stats, output_path + "/regions", prof);
    }
    prof.pop();
}

template<class M, class E, class V, class F, class P> inline
void SegmentedMesh::output(Polyhedralmesh<M,E,V,F,P> &m, Parameters &Par)
{
    prof.push("FESA::output");
    update_regions_map(m, polys_in_subregion);

    // print field and isovalues info
    print_field_csv(m, output_path + "/output_field.csv", prof);
    print_isovals_csv(m, percentiles, isovals, output_path + "/isovalues.txt", prof);

    // save mesh
    std::string domain_dir = output_path + "/domain";
    open_directory(domain_dir);
    m.save((domain_dir + "/domain.mesh").c_str());

    // print verts and cells data (label and field) in a csv file
    print_verts_csv(m, domain_dir + "/domain_verts_data.csv", prof);
    print_cells_csv(m, subregion_region_map, domain_dir + "/domain_cells_data.csv", prof);

    // print general statistics
    Statistics stats;
    stats.init(percentiles, isovals, prof);
    stats.compute_stats(m, polys_in_subregion);

    if (Par.get_OUT_LEVEL() > 1) {
        if (Par.get_OUT_LEVEL() > 2) {
            // save a mesh for each iso-subregion
            stats.init(percentiles, isovals, prof);
            print_regions(m, subregion_region_map, polys_in_subregion, stats, output_path + "/subregions", prof);
        }
        // save a mesh for each iso-region
        restore_original_labels(m);
        print_regions(m, subregion_region_map, polys_in_region, stats, output_path + "/regions", prof);
    }
    prof.pop();
}

/**********************************************************************/
/*************************** PRIVATE METHODS **************************/
/**********************************************************************/

template<class M, class V,  class E, class P> inline
void SegmentedMesh::load_field(AbstractMesh<M,V,E,P> &m, const std::string field_path, const int field_type)
{
    prof.push("Input field:  " + field_path);
    // load the field from file
    std::ifstream fp(field_path.c_str());
    assert(fp.is_open());
    double f;
    int i = 0;
    while (fp >> f) {
        field[i] = f + field_correction;
        ++i;
    }
    fp.close();
    assert(m.num_polys() <= field.size());

    if (field_type == 1) {
        assert(field.size() == m.num_polys());
        for (uint pid=0; pid<m.num_polys(); ++pid) {
            m.poly_data(pid).fvalue = field.at(pid);
        }
    } else if (field_type == 2) {
        assert(field.size() == m.num_verts());
        for (uint vid=0; vid<m.num_verts(); ++vid) {
            m.vert_data(vid).fvalue = field.at(vid);
        }
    }

    prof.pop();
}

template<class M, class V,  class E, class P> inline
void SegmentedMesh::extend_field_to_verts(AbstractMesh<M,V,E,P> &m)
{
    // set vert_data.uvw.u from the weighted average of poly.data.fvalue
    for (uint vid=0; vid<m.num_verts(); ++vid) {
        double value = 0., mass  = 0.;
        for (uint pid : m.adj_v2p(vid)) {
            double coeff = m.poly_mass(pid) / m.verts_per_poly(pid);
            value += field.at(pid) * coeff;
            mass  += coeff;
        }
        if (mass < 1e-4)
            m.vert_data(vid).fvalue = field.at(m.adj_v2p(vid).front());
        // assert(mass > 0);
        else
            m.vert_data(vid).fvalue = value / mass;
    }
}

template<class M, class V,  class E, class P> inline
void SegmentedMesh::extend_field_to_cells(AbstractMesh<M,V,E,P> &m)
{
    // set poly.data.fvalue from the average of vert_data.uvw.u
    for (uint pid=0; pid<m.num_polys(); ++pid) {
        double value = 0.;
        for (uint vid : m.adj_p2v(pid)) {
            value += m.vert_data(vid).fvalue;
        }
        m.poly_data(pid).fvalue = value / m.verts_per_poly(pid);
    }
}

void SegmentedMesh::load_global_field(const std::string field_path)
{
    prof.push("Global field: " + field_path);
    std::ifstream fp(field_path.c_str());
    assert(fp.is_open());
    double f;
    while (fp >> f) {
        global_field.push_back(f + field_correction);
    }
    fp.close();

    auto minmax = minmax_element(global_field.begin(), global_field.end());
    field_min = *minmax.first;
    field_max = *minmax.second;

    if (field_min < 0) {
        cout << "WARNING: the field has negative values! Field_min: " << field_min << endl << flush;
    }
    prof.pop();
}

template<class M, class V,  class E, class P> inline
void SegmentedMesh::color_mesh(AbstractMesh<M,V,E,P> &m) {
    for (uint pid=0; pid<m.num_polys(); ++pid) {
        // double c = (field.at(pid) - field_min) / (field_max - field_min); // input field
        double c = (m.poly_data(pid).fvalue - field_min) / (field_max - field_min); // mean field
        m.poly_data(pid).color = Color::red_white_blue_ramp_01(1. - c);
    }

    for (uint vid=0; vid<m.num_verts(); ++vid) {
        double c = (m.vert_data(vid).fvalue - field_min) / (field_max - field_min);
        m.vert_data(vid).color = Color::red_white_blue_ramp_01(1. - c);
        m.vert_data(vid).uvw[0] = m.vert_data(vid).fvalue; // for drawing isocontours
    }
}

/**********************************************************************/

void SegmentedMesh::compute_isovalues(const int n_regions, const int input_type, const std::vector<double> &input_vals)
{
    prof.push("Compute " + std::to_string(n_regions) + " iso-regions");

    switch (input_type) {
        case 1: { // equi-spaced isovals
            double step = abs(field_max - field_min) / n_regions;
            for (uint i = 0; i < n_regions; ++i) {
                double isoval = field_min + step * i;
                isovals.push_back(isoval);
                double perc = inv_percentile(global_field, isoval);
                percentiles.push_back(perc);
            }
            isovals.push_back(field_max);
            percentiles.push_back(100);
            break;
        }
        case 2: { // percentile isovals
            assert(input_vals.front()==0 && input_vals.back()==100);
            for (uint i = 0; i < input_vals.size()-1; ++i) {
                assert(input_vals.at(i) < input_vals.at(i+1) && "percentiles are not ordered");
            }
            percentiles = input_vals;
            isovals.push_back(field_min);
            for (int i=1; i<input_vals.size()-1; ++i) {
                int perc = input_vals.at(i);
                double isoval = percentile(global_field, perc);
                isovals.push_back(isoval);
            }
            isovals.push_back(field_max);
            break;
        }
        case 3: { // manually set isovals
            for (uint i = 0; i < input_vals.size()-1; ++i) {
                assert(input_vals.at(i) < input_vals.at(i+1) && "isovalues are not ordered");
            }
            for (double isoval : input_vals) {
                double perc = inv_percentile(global_field, isoval);
                percentiles.push_back(perc);
            }
            isovals = input_vals;
            break;
        }
        default: {
            assert(false && "unknown isovalues type");
            break;
        }
    }
    assert(isovals.front() == field_min && "first isoval != field_min");
    assert(isovals.back() == field_max && "last isoval != field_min");
    assert(isovals.size() == percentiles.size() && "#isovalues != #percentiles");
    assert((int)isovals.size() == n_regions + 1 && "#isovalues != #regions+1");
    prof.pop();
}

/**********************************************************************/

template<class M, class V,  class E, class P> inline
void SegmentedMesh::compute_isoregions(AbstractMesh<M,V,E,P> &m, const bool DENOISE)
{
    // assign the same label to all the polys in the same iso-region
    for (uint pid = 0; pid < m.num_polys(); ++pid) {
        bool found = false;
        double val = DENOISE ? poly_ring_mean(m, field, pid) : m.poly_data(pid).fvalue;
        assert(field_min <= val && val <= field_max);
        for (int lid = 0; lid < (int)isovals.size() - 1; ++lid) {
            if (isovals.at(lid) <= val && val < isovals.at(lid + 1)) {
                m.poly_data(pid).label = lid;
                found = true;
                break;
            }
        }
        if (val == isovals.back()) {
            m.poly_data(pid).label = isovals.size()-2;
            found = true;
        }
        assert(found);
    }
    update_regions_map(m, polys_in_region);
}

template<class M, class V,  class E, class P> inline
void SegmentedMesh::compute_isoregions2(AbstractMesh<M,V,E,P> &m, const bool DENOISE)
{
    // assign the same label to all the polys in the same iso-region
    for (uint pid = 0; pid < m.num_polys(); ++pid) {
        vector<int> vert_labels(isovals.size() - 1);
        for (uint vid : m.adj_p2v(pid)) {
            bool found = false;
            double val = DENOISE ? vert_ring_mean(m, field, vid) : m.vert_data(vid).fvalue;
            for (int lid = 0; lid < (int)isovals.size() - 1; ++lid) {
                if (isovals.at(lid) <= val && val < isovals.at(lid + 1)) {
                    ++vert_labels.at(lid);
                    found = true;
                    break;
                }
            }
            if (val == isovals.back()) {
                ++vert_labels.at(isovals.size()-2);
                found = true;
            }
            assert(found);
        }
        auto result = max_element(vert_labels.begin(), vert_labels.end());
        uint lid = distance(vert_labels.begin(), result);
        m.poly_data(pid).label = lid;
    }
    update_regions_map(m, polys_in_region);
}

/**********************************************************************/

template<class M, class V,  class E, class P> inline
void SegmentedMesh::separate_ccs(AbstractMesh<M,V,E,P> &m)
{
    prof.push("Separate conn. comp.");
    subregion_region_map.clear();
    int tmp_label = std::max_element(polys_in_region.begin(), polys_in_region.end(),
                    [](const auto& l, const auto& r) { return l.first < r.first; })->first + 1;
    std::vector<bool> visited(m.num_polys(), false);

    for (auto &l : polys_in_region) {
        // mask all the polys with a label different from *l.first* (if mask[p] = true, bfs cannot pass through it)
        std::vector<bool> mask(m.num_polys(), true);
        for (uint pid : l.second) {
            mask.at(pid) = false;
        }

        for (uint pid : l.second) {
            if (visited.at(pid)) continue;

            // expand from *pid* and conquer its connected component
            std::unordered_set<uint> component;
            bfs_on_dual(m, pid, mask, component);

            // assign a temporary label to the polys in the connected component
            for (uint pid : component) {
                m.poly_data(pid).label = tmp_label;
                visited.at(pid) = true;
            }

            // update the map between tmp labels and old labels
            auto it = subregion_region_map.find(l.first);
            int old_label = (it == subregion_region_map.end()) ? l.first : it->second;
            subregion_region_map[tmp_label] = old_label;
            ++tmp_label;
        }
    }
    update_regions_map(m, polys_in_subregion);
    prof.pop(true, "\tfound " + std::to_string(polys_in_subregion.size()) + " components");
}

/**********************************************************************/

template<class M, class V,  class E, class P> inline
void SegmentedMesh::remove_small_labels(AbstractMesh<M,V,E,P> &m, const double CLEAN_THRESH)
{
    if (verbose) prof.push("   Remove small conn. comp.");
    // find the labels with area smaller than threshold
    std::unordered_set<int> labels_to_remove;
    for (auto &l : polys_in_subregion) {
        bool CRITERION_1 = CRIT_mass(m, l.second, CLEAN_THRESH);
        if (CRITERION_1) {
            labels_to_remove.insert(l.first);
        }
    }

    // replace each poly in *labels_to_remove* with the dominant neighboring label around it
    std::unordered_map<uint, int> new_polys_label;
    uint count = 0;
    for (int l : labels_to_remove) {
        // count the occurrences of each label (except *l*) around elements with label *l*
        std::unordered_map<int, int> label_counter;
        for (uint pid : polys_in_subregion[l]) {
            for (uint nbr : m.adj_p2p(pid)) {
                int neigh_l = m.poly_data(nbr).label;
                if (DOES_NOT_CONTAIN(labels_to_remove, neigh_l)) {
                    ++label_counter[neigh_l];
                }
            }
        }
        if (label_counter.empty()) continue;

        // if there is a dominant label, map it to pid, otherwise keep the old label
        std::map<int, std::vector<int>> label_counter_ordered;
        for (auto &l : label_counter) {
            // order the map by value by switching keys and values
            label_counter_ordered[l.second].push_back(l.first);
        }
        // assert((--label_counter_ordered.end())->second.size() == 1);  // more than one candidate
        for (uint pid : polys_in_subregion[l]) {
            new_polys_label[pid] = (--label_counter_ordered.end())->second.front();
        }
        ++count;
    }

    // assign the new labels
    for (auto &l : new_polys_label) {
        m.poly_data(l.first).label = l.second;
    }
    update_regions_map(m, polys_in_subregion);
    if (verbose) prof.pop(true, "\tremoved " + std::to_string(count) + " components");
}

/**********************************************************************/

template<class M, class V,  class E, class P> inline
uint SegmentedMesh::smooth_boundaries(AbstractMesh<M,V,E,P> &m)
{
    if (verbose) prof.push("   Smoothing of the boundaries");

    uint count = 0;
    vector<bool> visited(m.num_polys(), false);
    for (auto &l : polys_in_subregion) {
        for (uint pid : l.second) {
            if (visited.at(pid)) continue;

            // count how many different labels there are around pid,
            // and how often each of them occurs
            std::unordered_map<int, uint> label_counter;
            for(uint nbr : m.adj_p2p(pid)) {
                int nbr_l = m.poly_data(nbr).label;
                if (nbr_l != l.first) {
                    ++label_counter[nbr_l];
                }
            }
            if (label_counter.empty()) {
                visited.at(pid) = true;
                continue;
            }

            // order the map by value by switching keys and values:
            // now the dominant label is at the end of the list
            std::map<uint, int> label_counter_ordered;
            for (auto &it : label_counter) {
                label_counter_ordered[it.second] = it.first;
            }
            uint max_occurence = (--label_counter_ordered.end())->first;
            int l_dominant = label_counter_ordered[max_occurence];

            // if the dominant label occurs on more than half of the neighbors of pid,
            // assign pid this label
            if (max_occurence > m.verts_per_poly(pid) - max_occurence) {
                m.poly_data(pid).label = l_dominant;
                ++count;
            }
            visited.at(pid) = true;
        }
    }
    update_regions_map(m, polys_in_subregion);

    if (verbose) prof.pop(true, "\tswitched " + std::to_string(count) + " labels");
    return count;
}

/**********************************************************************/

template<class M, class V,  class E, class P> inline
void SegmentedMesh::restore_original_labels(AbstractMesh<M,V,E,P> &m)
{
    prof.push("Reassign original labels: ");
    for (uint pid=0; pid<m.num_polys(); ++pid) {
        int tmp_label = m.poly_data(pid).label;
        auto it = subregion_region_map.find(tmp_label);
        // find the label corresponding to *tmp_label* in *subregion_region_map*
        // if *tmp_label* is not in the map, then it is one of the old labels and we maintain it
        int old_label = (it == subregion_region_map.end()) ? tmp_label : it->second;
        m.poly_data(pid).label = old_label;
    }
    prof.pop(true, "\tfound " + std::to_string(m.polys_n_unique_labels()) + " labels");
}
