#include "contoured_mesh.h"
#include "read_parameters.h"

void ContouredMesh::run(Polygonmesh<> &m, Parameters &Par)
{
    verbose = Par.get_VERBOSE();
    prof.push("\nInput mesh:   " + Par.get_MESH_PATH());
    m.load(Par.get_MESH_PATH().c_str());
    prof.pop();

    load_global_field(Par.get_FIELD2_PATH());
    load_field(m, Par.get_FIELD1_PATH());
    color_mesh(m);
    // load_second_field(m, Par.get_FIELD2_PATH());

    compute_iso(Par.get_N_REGIONS(), Par.get_ISOVALS());
    if (Par.get_CUT_MESH() && m.mesh_type() == TRIMESH) {
        // cut_mesh(m); // needs the Trimesh class instead of Polygonmesh
    }
    insert_iso2(m);

    separate_ccs(m);

    prof.push("Remove small conn. comp.");
    int i = 0, n_labels_tmp = 0;
    int n_labels_init = labels_polys_map.size();
    double FILTER_THRESH = m.mesh_area() * Par.get_FILTER_THRESH() / 100.;
    while (n_labels_tmp != (int)labels_polys_map.size() && i < Par.get_N_ITER()) {
        n_labels_tmp = labels_polys_map.size();
        remove_small_labels(m, FILTER_THRESH);
        ++i;
    }
    if (n_labels_tmp == n_labels_init)
        cout << "ContouredMesh Warning: did not remove any label. Filter_thresh is too high?";
    prof.pop(true, "\tfound " + std::to_string(labels_polys_map.size()) + " components");

    prof.push("Smoothing of the boundaries");
    mark_edges(m);
    int j = 0, count = 1;
    while (count > 0 && j < Par.get_N_ITER()) {
        count = smooth_boundaries(m);
        ++j;
    }
    prof.pop();

    mark_edges(m);

    Statistics stats(Par.get_MESH_PATH(), Par.get_ISOVALS(), isovals, labels_polys_map);
    stats.compute_stats(m);
    stats.print_stats(m);
    // stats.print_regions_csv(m);
    stats.print_regions_shp(m);

    restore_original_labels(m);

    if (Par.get_OUT_FORMAT() == "-")
        Par.set_OUT_FORMAT("obj");
    if (Par.get_OUT_FORMAT() != "0")
        save_mesh(m, Par.get_MESH_PATH(), Par.get_OUT_FORMAT());
}

void ContouredMesh::run(Polyhedralmesh<> &m, Parameters &Par)
{
    verbose = Par.get_VERBOSE();
    prof.push("Input mesh:   " + Par.get_MESH_PATH());
    m.load(Par.get_MESH_PATH().c_str());
    prof.pop();

    load_global_field(Par.get_FIELD2_PATH());
    load_field(m, Par.get_FIELD1_PATH());
    color_mesh(m);

    compute_iso(Par.get_N_REGIONS(), Par.get_ISOVALS());
    if (Par.get_CUT_MESH()) {
        // cut_mesh(m); // needs the Tetmesh class instead of Polyhedralmesh
    }
    insert_iso2(m);

    separate_ccs(m);

    prof.push("Remove small conn. comp.");
    int i = 0, n_labels_tmp = 0;
    int n_labels_init = labels_polys_map.size();
    double FILTER_THRESH = m.mesh_volume() * Par.get_FILTER_THRESH() / 100.;
    while (n_labels_tmp != (int)labels_polys_map.size() && i < Par.get_N_ITER()) {
        n_labels_tmp = labels_polys_map.size();
        remove_small_labels(m, FILTER_THRESH);
        ++i;
    }
    if (n_labels_tmp == n_labels_init)
        cout << "ContouredMesh Warning: did not remove any label. Filter_thresh is too high?";
    prof.pop(true, "\tfound " + std::to_string(labels_polys_map.size()) + " components");

    prof.push("Smoothing of the boundaries");
    mark_edges(m);
    int j = 0, count = 1;
    while (count > 0 && j < Par.get_N_ITER()) {
        count = smooth_boundaries(m);
        ++j;
    }
    prof.pop();

    mark_edges(m);

    Statistics stats(Par.get_MESH_PATH(), Par.get_ISOVALS(), isovals, labels_polys_map);
    stats.compute_stats(m);
    stats.print_stats(m);
    // stats.print_regions_csv(m);
    // stats.print_regions_shp(m);

    restore_original_labels(m);

    if (Par.get_OUT_FORMAT() == "-")
        Par.set_OUT_FORMAT("vtk");
    if (Par.get_OUT_FORMAT() != "0")
        save_mesh(m, Par.get_MESH_PATH(), Par.get_OUT_FORMAT());
}

/**********************************************************************/

void ContouredMesh::load_global_field(const std::string field_path)
{
    prof.push("Global field: " + field_path);
    std::ifstream fp(field_path.c_str());
    assert(fp.is_open());
    double f;
    while (fp >> f) {
        field_global.push_back(f + field_correction);
    }
    fp.close();

    auto minmax = minmax_element(field_global.begin(), field_global.end());
    field_min = *minmax.first;
    field_max = *minmax.second;

    if (field_min < 0) {
        cout << "WARNING: the field has negative values! Field_min: " << field_min << endl << flush;
    }
    prof.pop();
}

template<class M, class V, class E, class P> inline
void ContouredMesh::color_mesh(AbstractMesh<M,E,V,P> &m) {
    PARALLEL_FOR(0, m.num_polys(), 1000, [this, &m](int pid) {
        double c = (field.at(pid) - field_min) / (field_max - field_min);
        m.poly_data(pid).color = Color::red_white_blue_ramp_01(1. - c);
    });

    PARALLEL_FOR(0, m.num_verts(), 1000, [this, &m](int vid) {
        double c = (m.vert_data(vid).uvw.u() - field_min) / (field_max - field_min);
        m.vert_data(vid).color = Color::red_white_blue_ramp_01(1. - c);
    });
}

template<class M, class V, class E, class P> inline
void ContouredMesh::load_field(AbstractMesh<M,E,V,P> &m, const std::string field_path)
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
    PARALLEL_FOR(0, m.num_polys(), 1000, [this, &m](int pid) {
       m.poly_data(pid).quality = field.at(pid);
    });

    // set vert_data.uvw.u and vert_data.color from the weighted average of polys quality
    PARALLEL_FOR(0, m.num_verts(), 1000, [this, &m](int vid) {
        double value = 0., mass  = 0.;
        for (uint pid : m.adj_v2p(vid)) {
            double coeff = m.poly_mass(pid) / m.verts_per_poly(pid);
            value += field.at(pid) * coeff;
            mass  += coeff;
        }
        value /= mass;
        m.vert_data(vid).uvw.u() = value;
    });
    prof.pop();
}

template<class M, class V, class E, class P> inline
void ContouredMesh::load_second_field(AbstractMesh<M,E,V,P> &m, const std::string field_path)
{
    prof.push("Input field2: \t" + field_path);
    // load the field from file
    std::ifstream fp(field_path.c_str());
    assert(fp.is_open());
    double f;
    int _pid = 0;
    while (fp >> f) {
        field2[_pid] = f;
        ++_pid;
    }
    fp.close();
    assert(m.num_polys() <= field2.size());

    // store the field value in each poly
    std::vector<double> poly_heights(m.num_polys());
    PARALLEL_FOR(0, poly_heights.size(), 1000, [this, &poly_heights](int pid) {
        poly_heights.at(pid) = field2.at(pid);
    });

    // set vert.y to a weighted average of the neighboring poly values
    PARALLEL_FOR(0, m.num_verts(), 1000, [poly_heights, &m](int vid) {
        double value = 0., mass  = 0.;
        for (uint pid : m.adj_v2p(vid)) {
            double coeff = m.poly_mass(pid) / m.verts_per_poly(pid);
            value += poly_heights.at(pid) * coeff;
            mass  += coeff;
        }
        value /= mass;
        m.vert(vid).y() = value; //exp(value);
    });
    prof.pop();
}

/**********************************************************************/

// compute the *perc*-th percentile relative to the values in *data*
double percentile(const vector<double> &data, const double perc)
{
    assert(!data.empty());
    vector<double> sortedData = data;
    sort(sortedData.begin(), sortedData.end());
    int n = sortedData.size()-1;

    int index = perc * n / 100;
    if (index == n)
        return sortedData[index];

    double interpolation = perc * n / 100 - index;
    return sortedData[index] + interpolation * (sortedData[index+1] - sortedData[index]);
}

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
            double perc = percentile(field_global, val);
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

template<class M, class V, class E, class P> inline
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
    update_labels_polys_map(m);
}

// compute the mean value of *field* in the 1-ring of poly *pid*
template<class M, class V, class E, class P> inline
double poly_ring_mean(AbstractMesh<M,E,V,P> &m, std::map<uint,double> &field, const uint pid)
{
    vector<uint> ring;
    for (uint vid : m.adj_p2v(pid)) {
        vector<uint> nbr = m.adj_v2p(vid);
        ring.insert(ring.end(), nbr.begin(), nbr.end());
    }
    REMOVE_DUPLICATES_FROM_VEC(ring);
    vector<double> ring_values;
    for (uint pid : ring) {
        ring_values.push_back(field[pid]);
    }
    auto [min, max] = minmax_element(ring_values.begin(), ring_values.end());
    return sqrt(fabs(*min * *max));
}

template<class M, class V, class E, class P> inline
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
    update_labels_polys_map(m);
}

// compute the mean value of *field* in the 1-ring of vertex *vid*
template<class M, class V, class E, class P> inline
double vert_ring_mean(AbstractMesh<M,E,V,P> &m, std::map<uint,double> &field, const uint vid)
{
    vector<double> ring_values;
    for (uint pid : m.adj_v2p(vid)) {
        ring_values.push_back(field[pid]);
    }
    auto [min, max] = minmax_element(ring_values.begin(), ring_values.end());
    return sqrt(fabs(*min * *max));
}

template<class M, class V, class E, class P> inline
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
        m.poly_data(pid).label = distance(vert_labels.begin(), result);
    }
    update_labels_polys_map(m);
}

/**********************************************************************/

template<class M, class V, class E, class P> inline
void ContouredMesh::separate_ccs(AbstractMesh<M,E,V,P> &m)
{
    prof.push("Separate conn. comp.");
    tmp_labels_map.clear();
    int tmp_label = std::max_element(labels_polys_map.begin(), labels_polys_map.end(),
                    [](const auto& l, const auto& r) { return l.first < r.first; })->first + 1;
    std::vector<bool> visited(m.num_polys(), false);

    for (auto &l : labels_polys_map) {
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
            auto it = tmp_labels_map.find(l.first);
            int old_label = (it == tmp_labels_map.end()) ? l.first : it->second;
            tmp_labels_map[tmp_label] = old_label;
            ++tmp_label;
        }
    }
    update_labels_polys_map(m);
    prof.pop(true, "\tfound " + std::to_string(labels_polys_map.size()) + " components");
}

/**********************************************************************/

template<class M, class V, class E, class P> inline
void ContouredMesh::remove_small_labels(AbstractMesh<M,E,V,P> &m, const double FILTER_THRESH)
{
    if (verbose) prof.push("   Remove small conn. comp.");
    // find the labels with area smaller than threshold
    std::unordered_set<int> labels_to_remove;
    for (auto &l : labels_polys_map) {
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
        for (uint pid : labels_polys_map[l]) {
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
        for (uint pid : labels_polys_map[l]) {
            new_polys_label[pid] = (--label_counter_ordered.end())->second.front();
        }
        ++count;
    }

    // assign the new labels
    for (auto &l : new_polys_label) {
        m.poly_data(l.first).label = l.second;
    }
    update_labels_polys_map(m);
    if (verbose) prof.pop(true, "\tremoved " + std::to_string(count) + " components");
}

/**********************************************************************/

template<class M, class V, class E, class P> inline
uint ContouredMesh::smooth_boundaries(AbstractMesh<M,E,V,P> &m)
{
    if (verbose) prof.push("   Smoothing of the boundaries");

    uint count = 0;
    vector<bool> visited(m.num_polys(), false);
    for (auto &l : labels_polys_map) {
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
    update_labels_polys_map(m);

    if (verbose) prof.pop(true, "\tswitched " + std::to_string(count) + " labels");
    return count;
}

/**********************************************************************/

template<class M, class V, class E, class P> inline
void ContouredMesh::restore_original_labels(AbstractMesh<M,E,V,P> &m)
{
    prof.push("Reassign original labels: ");
    PARALLEL_FOR(0, m.num_polys(), 1000, [this, &m](int pid) {
        int tmp_label = m.poly_data(pid).label;
        auto it = tmp_labels_map.find(tmp_label);
        // find the label corresponding to *tmp_label* in *tmp_labels_map*
        // if *tmp_label* is not in the map, then it is one of the old labels and we maintain it
        int old_label = (it == tmp_labels_map.end()) ? tmp_label : it->second;
        m.poly_data(pid).label = old_label;
    });
    prof.pop(true, "\tfound " + std::to_string(m.polys_n_unique_labels()) + " labels");
}

/**********************************************************************/

template<class M, class V, class E, class P> inline
void ContouredMesh::mark_edges(AbstractMesh<M,E,V,P> &m)
{
    PARALLEL_FOR(0, m.num_edges(), 1000, [&m](int eid) {
        m.edge_data(eid).flags[MARKED] = false;
        std::vector<uint> polys = m.adj_e2p(eid);
        if (m.poly_data(polys.front()).label != m.poly_data(polys.back()).label) {
            m.edge_data(eid).flags[MARKED] = true;
        }
    });
}

void ContouredMesh::mark_faces(Tetmesh<> &m)
{
    PARALLEL_FOR(0, m.num_faces(), 1000, [&m](int fid) {
        m.face_data(fid).flags[MARKED] = false;
        for (uint eid : m.adj_f2e(fid)) {
            m.edge_data(eid).flags[MARKED] = false;
        }

        std::vector<uint> polys = m.adj_f2p(fid);
        if (m.poly_data(polys.front()).label != m.poly_data(polys.back()).label) {
            m.face_data(fid).flags[MARKED] = true;
            for (uint eid : m.adj_f2e(fid)) {
                m.edge_data(eid).flags[MARKED] = true;
            }
        }
    });
}

/**********************************************************************/

template<class M, class V, class E, class P> inline
void ContouredMesh::update_labels_polys_map(AbstractMesh<M,E,V,P> &m)
{
    labels_polys_map.clear();
    for (uint pid = 0; pid < m.num_polys(); ++pid) {
        labels_polys_map[m.poly_data(pid).label].push_back(pid);
    }
}

/**********************************************************************/

template<class M, class V, class E, class P> inline
void ContouredMesh::save_mesh(AbstractMesh<M,E,V,P> &m,
                              const std::string output_path, const std::string MESH_FORMAT)
{
    prof.push("Save mesh");
    // save the mesh
    size_t lastindex = std::string(output_path).find_last_of(".");
    std::string mesh_path = std::string(output_path).substr(0, lastindex) + "_out." + MESH_FORMAT;
    // m.save(mesh_path.c_str());

    // save the cells labels
    std::string labels_path = std::string(output_path).substr(0, lastindex) + "_labels.csv";
    std::ofstream fp;
    fp.open(labels_path.c_str());
    assert(fp.is_open());
    fp << "SCALAR_FIELD " << m.num_polys() << "\n";
    for (uint pid = 0; pid < m.num_polys(); ++pid) {
        fp << m.poly_data(pid).label << std::endl;
    }
    fp.close();
    prof.pop();
}
