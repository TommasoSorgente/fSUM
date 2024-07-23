#ifndef AUXILIARY_H
#define AUXILIARY_H
 
#include <filesystem>
#include <shapefil.h>

#include <cinolib/meshes/meshes.h>

using namespace cinolib;
namespace fs = std::filesystem;

/* Auxiliary functions */

/**********************************************************************/

// compute the *perc*-th percentile relative to the values in *data*
double percentile(const std::vector<double> &data, const double perc)
{
    assert(!data.empty());
    std::vector<double> sortedData = data;
    sort(sortedData.begin(), sortedData.end());
    int n = sortedData.size()-1;

    int index = perc * n / 100;
    if (index == n)
        return sortedData[index];

    double interpolation = perc * n / 100 - index;
    return sortedData[index] + interpolation * (sortedData[index+1] - sortedData[index]);
}

/**********************************************************************/

// compute the mean value of *field* in the 1-ring of poly *pid*
template<class M, class V, class E, class P> inline
double poly_ring_mean(AbstractMesh<M,E,V,P> &m, std::map<uint,double> &field, const uint pid)
{
    std::vector<uint> ring;
    for (uint vid : m.adj_p2v(pid)) {
        std::vector<uint> nbr = m.adj_v2p(vid);
        ring.insert(ring.end(), nbr.begin(), nbr.end());
    }
    REMOVE_DUPLICATES_FROM_VEC(ring);
    std::vector<double> ring_values;
    for (uint pid : ring) {
        ring_values.push_back(field[pid]);
    }
    auto [min, max] = minmax_element(ring_values.begin(), ring_values.end());
    return sqrt(fabs(*min * *max));
}

// compute the mean value of *field* in the 1-ring of vertex *vid*
template<class M, class V, class E, class P> inline
double vert_ring_mean(AbstractMesh<M,E,V,P> &m, std::map<uint,double> &field, const uint vid)
{
    std::vector<double> ring_values;
    for (uint pid : m.adj_v2p(vid)) {
        ring_values.push_back(field[pid]);
    }
    auto [min, max] = minmax_element(ring_values.begin(), ring_values.end());
    return sqrt(fabs(*min * *max));
}

/**********************************************************************/

/* mark edges or faces among cells with different labels */
template<class M, class V, class E, class P> inline
void mark_edges(AbstractMesh<M,E,V,P> &m)
{
    PARALLEL_FOR(0, m.num_edges(), 1000, [&m](int eid) {
        m.edge_data(eid).flags[MARKED] = false;
        std::vector<uint> polys = m.adj_e2p(eid);
        if (m.poly_data(polys.front()).label != m.poly_data(polys.back()).label) {
            m.edge_data(eid).flags[MARKED] = true;
        }
    });
}

void mark_faces(Tetmesh<> &m)
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

void open_directory(const std::string &path) {
    if (!fs::exists(path)) {
        // if the folder does not exist, create it
        try {
            fs::create_directories(path);
        } catch (const std::exception &e) {
            std::cerr << "Error creating folder: " << e.what() << std::endl;
            assert(false);
        }
    } else {
        // if the folder already exists, delete its content
        for (const auto &entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file()) {
                fs::remove(entry.path());
            } else if (entry.is_directory()) {
                fs::remove_all(entry.path());
            }
        }
    }
}

/**********************************************************************/

void extract_region(Polygonmesh<> &m, const uint label,
                    std::vector<std::vector<vec3d>> &verts) {
    Polygonmesh<> sub_m;
    export_cluster(m, label, sub_m);
    assert(connected_components(sub_m) == 1);

    // count how many times each boundary vertex should be visited
    std::vector<int> n_visits(sub_m.num_verts(), 0);
    uint total_visits = 0;
    for (uint vid = 0; vid < sub_m.num_verts(); ++vid) {
        for (uint eid : sub_m.adj_v2e(vid)) {
            if (sub_m.edge_is_boundary(eid)) {
                ++n_visits.at(vid);
            }
        }
        assert(n_visits.at(vid) % 2 == 0);
        n_visits.at(vid) /= 2;
        total_visits += n_visits.at(vid);
    }

    std::vector<uint> boundary_verts = sub_m.get_boundary_vertices();
    int v_curr = boundary_verts.front(), v_prev = -1, v_next = -1;
    uint visits_count = 0;
    while (visits_count < total_visits) {
        // vertices in a connected component of the region boundary
        std::vector<vec3d> verts_cc;
        while (n_visits.at(v_curr) > 0) {
            // add v_curr to verts_cc
            verts_cc.push_back(sub_m.vert(v_curr));
            --n_visits.at(v_curr);

            // find v_next
            if (sub_m.vert_is_manifold(v_curr)) {
                // pick v_next on the boundary edge around v_curr which does not contain v_prev
                for (uint eid : sub_m.adj_v2e(v_curr)) {
                    if (sub_m.edge_is_boundary(eid) &&
                        !sub_m.edge_contains_vert(eid, v_prev)) {
                        v_next = sub_m.vert_opposite_to(eid, v_curr);
                        break;
                    }
                }
            } else {
                // if v_curr is not manifold, find v_next through adjacencies
                if (v_prev == -1) {
                    // if the non-manifold vert is at the beginning of the list, v_prev is not defined
                    // set v_prev as a random boundary vert around v_curr
                    for (uint eid : sub_m.adj_v2e(v_curr)) {
                        if (sub_m.edge_is_boundary(eid)) {
                            v_prev = sub_m.vert_opposite_to(eid, v_curr);
                            break;
                        }
                    }
                }
                uint eid = sub_m.edge_id(v_prev, v_curr);
                int pid = sub_m.adj_e2p(eid).front();
                while (pid != -1) {
                    for (uint ejd : sub_m.adj_p2e(pid)) {
                        if (ejd != eid && sub_m.edge_contains_vert(ejd, v_curr)) {
                            eid = ejd;
                            break;
                        }
                    }
                    pid = sub_m.poly_opposite_to(eid, pid);
                }
                assert(sub_m.edge_is_boundary(eid));
                v_next = sub_m.vert_opposite_to(eid, v_curr);
            }
            assert(v_next != v_curr && "extract_region: could not find v_next");

            // update v_prev and v_curr
            v_prev = v_curr;
            v_curr = v_next;
            ++visits_count;
        }
        assert(sub_m.vert(v_next) == verts_cc.front() && "extract_region: open sequence");
        verts.push_back(verts_cc);

        // if there are still boundary vertices to visit, jump to an unvisited one
        if (visits_count < total_visits) {
            for (uint vid : boundary_verts) {
                if (n_visits.at(vid) > 0) {
                    v_curr = vid;
                    v_prev = -1;
                    v_next = -1;
                    break;
                }
            }
            assert(v_next == -1 && "extract_region: failed to jump to another cc");
        }
    }
}

void print_regions_csv(Polygonmesh<> &m, const std::unordered_map<int, std::vector<uint>> &labels_polys_map, const std::vector<int> &label_vec, Profiler Prof) {
    Prof.push("Print Regions");
    std::string base_path = std::string(HOME_PATH) + "out/contours_csv";
    open_directory(base_path);
    //
    for (uint i = 0; i < labels_polys_map.size(); ++i) {
        uint pid0 = labels_polys_map.at(label_vec[i]).front();
        int label = m.poly_data(pid0).label;
        std::vector<std::vector<vec3d>> verts;
        extract_region(m, label, verts);

        std::string path = base_path + "/region_" + std::to_string(i) + ".csv";
        std::ofstream fp;
        fp.open(path.c_str());
        assert(fp.is_open());
        fp << std::fixed;
        fp.precision(2);

        // header
        fp << "# n_boundary_components, n_vertices_per_component, vertices" << std::endl;

        // number of connected components of the region boundary
        int n_comp = verts.size();
        fp << n_comp << std::endl;

        // number of vertices in each connected component
        for (uint j = 0; j < verts.size(); ++j) {
            fp << verts.at(j).size() << std::endl;
        }

        // list of vertices in each connected component
        for (std::vector<vec3d> &verts_cc : verts) {
            for (vec3d &v : verts_cc) {
                fp << v.x() << " " << v.y() << " " << v.z() << std::endl;
            }
        }
        fp.close();
    }
    Prof.pop();
}

/**********************************************************************/

void print_regions_shp(Polygonmesh<> &m, const std::unordered_map<int, std::vector<uint>> &labels_polys_map, const std::vector<int> &label_vec, const std::string base_path, Profiler Prof) {
    Prof.push("Print Regions Shapefile");
    std::string filename = base_path + "/contours.shp";
    SHPHandle hSHP = SHPCreate(filename.c_str(), SHPT_POLYGON);

    for (uint i = 0; i < labels_polys_map.size(); ++i) {
        uint pid0 = labels_polys_map.at(label_vec[i]).front();
        int label = m.poly_data(pid0).label;
        std::vector<std::vector<vec3d>> verts;
        extract_region(m, label, verts);

        for (uint j = 0; j < verts.size(); ++j) {
            uint n = verts.at(j).size();
            double x[n], y[n];
            uint count = 0;

            // polygon vertices coordinates
            for (vec3d &v : verts.at(j)) {
                x[count] = v.x();
                y[count] = v.y();
                ++count;
            }

            // create polygon
            SHPObject* obj = SHPCreateSimpleObject(SHPT_POLYGON, n, x, y, nullptr);

            // write polygon on file
            SHPWriteObject(hSHP, -1, obj);
            SHPDestroyObject(obj);
        }
    }
    SHPClose(hSHP);
    Prof.pop();
}

/**********************************************************************/

/* save the final mesh and export the cell labels in a csv file */
template<class M, class V, class E, class P> inline
void save_mesh(AbstractMesh<M,E,V,P> &m,
               const std::string output_path, const std::string MESH_FORMAT)
{
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
}

#endif // AUXILIARY_H
