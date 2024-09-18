#include "contoured_mesh.h"

/**********************************************************************/
/*************************** PUBLIC METHODS ***************************/
/**********************************************************************/

template<class M, class E, class V, class P> inline
void ContouredMesh::init(AbstractMesh<M,E,V,P> &m, Parameters &Par)
{
    prof.push("FESA::init");
    verbose = Par.get_VERBOSE();

    // load mesh and scalar fields
    m.load(Par.get_MESH_PATH().c_str());

    load_field(m, Par.get_FIELD1_PATH(), field);
    if (Par.get_FIELDG_PATH() != "") {
        load_global_field(Par.get_FIELDG_PATH());
    } else {
        for (uint i=0; i<field.size(); ++i) {
            global_field.push_back(field[i]);
        }
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
    print_field_csv(m, output_path + "/input_field.csv", prof);
    prof.pop();
}

template<class M, class E, class V, class P> inline
void ContouredMesh::segment(AbstractMesh<M,E,V,P> &m, Parameters &Par)
{
    prof.push("FESA::segment");
    compute_iso(Par.get_N_REGIONS(), Par.get_PERCENTILES());
    if (Par.get_CUT_MESH() && (m.mesh_type() == TRIMESH || m.mesh_type() == TETMESH)) {
        // cut_mesh(m); // needs the Tri/Tetmesh class instead of Polygon/Polyhedralmesh
    }

    insert_iso2(m);

    separate_ccs(m); // from now on, we switch from regions to subregions!

    prof.pop();
}

/**********************************************************************/

template<class M, class E, class V, class P> inline
void ContouredMesh::filter(AbstractMesh<M,E,V,P> &m, Parameters &Par)
{
    prof.push("FESA::filter");
    int i = 0, n_labels_tmp = 0;
    int n_labels_init = polys_in_subregion.size();

    // compute mesh total mass (area or volume)
    double mass = 0;
    for (uint pid=0; pid<m.num_polys(); ++pid) {
        mass += m.poly_mass(pid);
    }
    double FILTER_THRESH = mass * Par.get_FILTER_THRESH() / 100.;

    while (n_labels_tmp != (int)polys_in_subregion.size() && i < Par.get_N_ITER()) {
        n_labels_tmp = polys_in_subregion.size();
        remove_small_labels(m, FILTER_THRESH);
        ++i;
    }

    if (n_labels_tmp == n_labels_init) {
        cout << "FESA Warning: did not remove any label. Filter_thresh is too high?";
    }
    prof.pop(true, "\tfound " + std::to_string(polys_in_subregion.size()) + " components");
}

/**********************************************************************/

template<class M, class E, class V, class P> inline
void ContouredMesh::smooth(AbstractMesh<M,E,V,P> &m, Parameters &Par)
{
    prof.push("FESA::smooth");
    int j = 0, count = 1;
    while (count > 0 && j < Par.get_N_ITER()) {
        count = smooth_boundaries(m);
        ++j;
    }
    prof.pop();
}

/**********************************************************************/

void ContouredMesh::output(Polygonmesh<> &m, Parameters &Par)
{
    prof.push("FESA::output");
    update_regions_map(m, polys_in_subregion);

    // print field and isovalues info
    print_field_csv(m, output_path + "/output_field.csv", prof);
    print_isovals_csv(m, Par.get_PERCENTILES(), isovals, output_path + "/isovalues.txt", prof);

    // save mesh
    std::string domain_dir = output_path + "/domain";
    open_directory(domain_dir);
    m.save((domain_dir + "/domain.obj").c_str());

    // print verts and cells data (label and field) in a csv file
    print_verts_csv(m, domain_dir + "/domain_verts_data.csv", prof);
    print_cells_csv(m, subregion_region_map, domain_dir + "/domain_cells_data.csv", prof);

    // print general statistics
    Statistics stats;
    stats.init(Par.get_PERCENTILES(), isovals, prof);
    stats.compute_stats(m, polys_in_subregion);
    stats.print_global_stats(m, Par.get_FILTER_THRESH(), output_path + "/global_stats.txt");
    print_misclassifications(m, isovals, subregion_region_map, output_path + "/misclassifications.csv", prof);

    // generate a shapefile containing the iso-regions boundaries
    print_regions_shp(m, stats, polys_in_subregion, domain_dir + "/domain", prof);

    // save a mesh for each iso-subregion
    stats.init(Par.get_PERCENTILES(), isovals, prof);
    print_regions(m, subregion_region_map, polys_in_subregion, stats, output_path + "/subregions", prof);

    // save a mesh for each iso-region
    restore_original_labels(m);
    print_regions(m, subregion_region_map, polys_in_region, stats, output_path + "/regions", prof);

    prof.pop();
}

void ContouredMesh::output(Polyhedralmesh<> &m, Parameters &Par)
{
    prof.push("FESA::output");
    update_regions_map(m, polys_in_subregion);

    // print field and isovalues info
    print_field_csv(m, output_path + "/output_field.csv", prof);
    print_isovals_csv(m, Par.get_PERCENTILES(), isovals, output_path + "/isovalues.txt", prof);

    // save mesh
    std::string domain_dir = output_path + "/domain";
    open_directory(domain_dir);
    m.save((domain_dir + "/domain.vtk").c_str());

    // print verts and cells data (label and field) in a csv file
    print_verts_csv(m, domain_dir + "/domain_verts_data.csv", prof);
    print_cells_csv(m, subregion_region_map, domain_dir + "/domain_cells_data.csv", prof);

    // print general statistics
    Statistics stats;
    stats.init(Par.get_PERCENTILES(), isovals, prof);
    stats.compute_stats(m, polys_in_subregion);
    stats.print_global_stats(m, Par.get_FILTER_THRESH(), output_path + "/global_stats.txt");
    print_misclassifications(m, isovals, subregion_region_map, output_path + "/misclassifications.csv", prof);

    // save a mesh for each iso-subregion
    stats.init(Par.get_PERCENTILES(), isovals, prof);
    print_regions(m, subregion_region_map, polys_in_subregion, stats, output_path + "/subregions", prof);

    // save a mesh for each iso-region
    restore_original_labels(m);
    print_regions(m, subregion_region_map, polys_in_region, stats, output_path + "/regions", prof);

    prof.pop();
}

/**********************************************************************/
/*************************** PRIVATE METHODS **************************/
/**********************************************************************/

template<class M, class E, class V, class P> inline
void ContouredMesh::load_field(AbstractMesh<M,E,V,P> &m, const std::string field_path, std::map<uint,double> &field)
{
    prof.push("Input field:  " + field_path);
    // load the field from file
    std::ifstream fp(field_path.c_str());
    assert(fp.is_open());
    double f;
    int _pid = 0;
    while (fp >> f) {
        field[_pid] = f + field_correction;
        ++_pid;
    }
    fp.close();
    assert(m.num_polys() <= field.size());

    // set poly_data.quality from the field values
    for (uint pid=0; pid<m.num_polys(); ++pid) {
        m.poly_data(pid).quality = field.at(pid);
    }

    // set vert_data.uvw.u and vert_data.color from the weighted average of polys quality
    for (uint vid=0; vid<m.num_verts(); ++vid) {
        double value = 0., mass  = 0.;
        for (uint pid : m.adj_v2p(vid)) {
            double coeff = m.poly_mass(pid) / m.verts_per_poly(pid);
            value += field.at(pid) * coeff;
            mass  += coeff;
        }
        value /= mass;
        m.vert_data(vid).uvw.u() = value;
    }
    prof.pop();
}

template<class M, class E, class V, class P> inline
void ContouredMesh::load_second_field(AbstractMesh<M,E,V,P> &m, const std::string field_path)
{
    prof.push("Input second_field: \t" + field_path);
    // load the field from file
    std::ifstream fp(field_path.c_str());
    assert(fp.is_open());
    double f;
    int _pid = 0;
    while (fp >> f) {
        second_field[_pid] = f;
        ++_pid;
    }
    fp.close();
    assert(m.num_polys() <= second_field.size());

    // store the field value in each poly
    std::vector<double> poly_heights(m.num_polys());
    for (uint pid=0; pid<m.num_polys(); ++pid) {
        poly_heights.at(pid) = second_field.at(pid);
    }

    // set vert.y to a weighted average of the neighboring poly values
    for (uint vid=0; vid<m.num_verts(); ++vid) {
        double value = 0., mass  = 0.;
        for (uint pid : m.adj_v2p(vid)) {
            double coeff = m.poly_mass(pid) / m.verts_per_poly(pid);
            value += poly_heights.at(pid) * coeff;
            mass  += coeff;
        }
        value /= mass;
        m.vert(vid).y() = value; //exp(value);
    }
    prof.pop();
}

void ContouredMesh::load_global_field(const std::string field_path)
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

template<class M, class E, class V, class P> inline
void ContouredMesh::color_mesh(AbstractMesh<M,E,V,P> &m) {
    for (uint pid=0; pid<m.num_polys(); ++pid) {
        double c = (field.at(pid) - field_min) / (field_max - field_min);
        m.poly_data(pid).color = Color::red_white_blue_ramp_01(1. - c);
    }

    for (uint vid=0; vid<m.num_verts(); ++vid) {
        double c = (m.vert_data(vid).uvw.u() - field_min) / (field_max - field_min);
        m.vert_data(vid).color = Color::red_white_blue_ramp_01(1. - c);
    }
}

/**********************************************************************/

void ContouredMesh::compute_iso(const int n_regions, const std::vector<double> &input_vals)
{
    prof.push("Compute " + std::to_string(n_regions) + " iso-regions");

    if (input_vals.empty()) { // equi-spaced isovals
        double step = abs(field_max - field_min) / n_regions;
        for (int lid = 0; lid < n_regions; ++lid) {
            isovals.push_back(field_min + step * lid);
        }
        isovals.push_back(field_max);
    }

    else { // percentile isovals
        assert(input_vals.front()==0 && input_vals.back()==100);
        for (uint i = 0; i < input_vals.size()-1; ++i) {
            assert(input_vals.at(i) < input_vals.at(i+1) && "input_vals are not ordered");
        }
        for (int val : input_vals) {
            double perc = percentile(global_field, val);
            isovals.push_back(perc);
            if (verbose) cout << "\tiso-value at percentile " << val << ": \t" << perc << endl;
        }
    }

    // else { // manually set isovals
    //     for (uint i = 0; i < input_vals.size()-1; ++i) {
    //         assert(input_vals.at(i) < input_vals.at(i+1) && "isovalues are not ordered");
    //     }
    //     isovals = input_vals;
    // }

    assert((int)isovals.size() == n_regions + 1 && "#isovalues != #regions+1");
    prof.pop();
}

/**********************************************************************/

void ContouredMesh::cut_mesh(Trimesh<> &m)
{
    // insert the isocontours in the mesh
    for (double v_value : isovals) {
        Isocontour<> iso(m, v_value);
        std::vector<uint> new_vids = iso.tessellate(m);
        for (uint vid : new_vids) {
            m.vert_data(vid).color = Color::red_white_blue_ramp_01(1. - v_value);
        }
    }
}

void ContouredMesh::cut_mesh(Tetmesh<> &m)
{
    // insert the isocontours in the mesh
    for (double value : isovals) {
        Isosurface<> iso(m, value);
        std::vector<uint> new_vids = iso.tessellate(m);
        for (uint vid : new_vids) {
            m.vert_data(vid).color = Color::red_white_blue_ramp_01(1. - value);
        }
    }
}

/**********************************************************************/

template<class M, class E, class V, class P> inline
void ContouredMesh::insert_iso1(AbstractMesh<M,E,V,P> &m)
{
    // set poly_data.quality from the vertices values
    for (uint pid = 0; pid < m.num_polys(); ++pid) {
        double p_value = 0.;
        for (uint vid : m.adj_p2v(pid)) {
            p_value += m.vert_data(vid).uvw.u();
        }
        m.poly_data(pid).quality = p_value / m.adj_p2v(pid).size();
    }
    // assign the same label to all the polys in the same iso-region
    for (uint pid = 0; pid < m.num_polys(); ++pid) {
        double val = m.poly_data(pid).quality;
        for (int lid = 0; lid < (int)isovals.size() - 1; ++lid) {
            if (isovals.at(lid) <= val && val < isovals.at(lid + 1)) {
                m.poly_data(pid).label = lid;
                break;
            }
            if (val == isovals.back()) {
                m.poly_data(pid).label = lid;
                break;
            }
        }
    }
    update_regions_map(m, polys_in_region);
}

template<class M, class E, class V, class P> inline
void ContouredMesh::insert_iso2(AbstractMesh<M,E,V,P> &m)
{
    // set poly_data.quality from the poly 1-ring
    for (uint pid = 0; pid < m.num_polys(); ++pid) {
       m.poly_data(pid).quality = poly_ring_mean(m, field, pid);
    }
    // assign the same label to all the polys in the same iso-region
    for (uint pid = 0; pid < m.num_polys(); ++pid) {
        bool found = false;
        double val = m.poly_data(pid).quality;
        for (int lid = 0; lid < (int)isovals.size() - 1; ++lid) {
            if (isovals.at(lid) <= val && val < isovals.at(lid + 1)) {
                m.poly_data(pid).label = lid;
                found = true;
                break;
            }
            if (val == isovals.back()) {
                m.poly_data(pid).label = lid;
                found = true;
                break;
            }
        }
        assert(found);
    }
    update_regions_map(m, polys_in_region);
}

template<class M, class E, class V, class P> inline
void ContouredMesh::insert_iso3(AbstractMesh<M,E,V,P> &m)
{
    // set vert_data.quality from the vert 1-ring
    for (uint vid = 0; vid < m.num_verts(); ++vid) {
        m.vert_data(vid).quality = vert_ring_mean(m, field, vid);
    }

    // assign the same label to all the polys in the same iso-region
    for (uint pid = 0; pid < m.num_polys(); ++pid) {
        vector<int> vert_labels(isovals.size() - 1);
        for (uint vid : m.adj_p2v(pid)) {
            bool found = false;
            double val = m.vert_data(vid).quality;
            for (int lid = 0; lid < (int)isovals.size() - 1; ++lid) {
                if (isovals.at(lid) <= val && val < isovals.at(lid + 1)) {
                    ++vert_labels.at(lid);
                    found = true;
                    break;
                }
                if (val == isovals.back()) {
                    ++vert_labels.at(lid);
                    found = true;
                    break;
                }
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

template<class M, class E, class V, class P> inline
void ContouredMesh::separate_ccs(AbstractMesh<M,E,V,P> &m)
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

template<class M, class E, class V, class P> inline
void ContouredMesh::remove_small_labels(AbstractMesh<M,E,V,P> &m, const double FILTER_THRESH)
{
    if (verbose) prof.push("   Remove small conn. comp.");
    // find the labels with area smaller than threshold
    std::unordered_set<int> labels_to_remove;
    for (auto &l : polys_in_subregion) {
        bool CRITERION_1 = CRIT_mass(m, l.second, FILTER_THRESH);
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

template<class M, class E, class V, class P> inline
uint ContouredMesh::smooth_boundaries(AbstractMesh<M,E,V,P> &m)
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

template<class M, class E, class V, class P> inline
void ContouredMesh::restore_original_labels(AbstractMesh<M,E,V,P> &m)
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
