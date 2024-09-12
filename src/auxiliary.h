#ifndef AUXILIARY_H
#define AUXILIARY_H
 
#include <filesystem>
#include <shapefil.h>

#include <cinolib/meshes/meshes.h>
#include <cinolib/export_cluster.h>
#include <cinolib/profiler.h>
#include <cinolib/connected_components.h>

#include "statistics.h"

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

/* create a map <label, all polys with that label> */
template<class M, class V, class E, class P> inline
void update_labels_polys_map(AbstractMesh<M,E,V,P> &m, std::unordered_map<int, std::vector<uint>> &labels_polys_map)
{
    labels_polys_map.clear();
    for (uint pid = 0; pid < m.num_polys(); ++pid) {
        labels_polys_map[m.poly_data(pid).label].push_back(pid);
    }
}


/**********************************************************************/

/* create the directory if it does not exist, otherwise delete its content */
void open_directory(const std::string &path, bool erase = true) {
    if (!fs::exists(path)) {
        // if the folder does not exist, create it
        try {
            fs::create_directories(path);
        } catch (const std::exception &e) {
            std::cerr << "Error creating folder: " << e.what() << std::endl;
            assert(false);
        }
    } else if (erase) {
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

/* extract the ordered sequence of vertices along the boundary of the region with label *label* */
void extract_region(Polygonmesh<> &m, const uint label,
                    std::vector<std::vector<uint>> &verts) {
    Polygonmesh<> sub_m;
    std::unordered_map<uint,uint> m2subm_vmap;
    std::unordered_map<uint,uint> subm2m_vmap;
    export_cluster(m, label, sub_m, m2subm_vmap, subm2m_vmap);
    // assert(connected_components(sub_m) == 1);

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
        std::vector<uint> verts_cc;
        while (n_visits.at(v_curr) > 0) {
            // add v_curr to verts_cc
            verts_cc.push_back(subm2m_vmap[v_curr]);
            // verts_cc.push_back(sub_m.vert(v_curr));
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
        // assert(sub_m.vert(v_next) == verts_cc.front() && "extract_region: open sequence");
        // assert(subm2m_vmap[v_next] == verts_cc.front() && "extract_region: open sequence");
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

/**********************************************************************/

/* print isovalues info */
template<class M, class V, class E, class P> inline
void print_isovals_csv(AbstractMesh<M,E,V,P> &m, const std::vector<double> &percentiles, const std::vector<double> &isovals, const std::string filename, Profiler Prof) {
    Prof.push("Print Isovalues");
    std::ofstream fp;
    fp.open(filename.c_str());
    assert(fp.is_open());
    fp << std::fixed;
    fp.precision(2);

    fp << "# percentile, isovalue, label" << std::endl;
    for (uint i=0; i<percentiles.size(); ++i) {
        fp << percentiles.at(i) << ", " << isovals.at(i);
        if (i < isovals.size() - 1)
            fp << ", " << i;
        fp << std::endl;
    }
    fp.close();
    Prof.pop();
}

/**********************************************************************/

/* print the vertices of the regions on separate csv files */
template<class M, class V, class E, class P> inline
void print_verts_csv(AbstractMesh<M,E,V,P> &m, const std::string filename, Profiler Prof) {
    Prof.push("Print Region Vertices");
    std::ofstream fp;
    fp.open(filename.c_str());
    assert(fp.is_open());
    fp << std::fixed;
    fp.precision(2);

    fp << "# vertex id, x, y, z, field" << std::endl;
    for (uint vid=0; vid<m.num_verts(); ++vid) {
        vec3d v = m.vert(vid);
        fp << vid << ", " << v.x() << ", " << v.y() << ", " << v.z() << ", " << m.vert_data(vid).uvw.u() << std::endl;
    }
    fp.close();
    Prof.pop();
}

/**********************************************************************/

/* print the cells of the regions on separate csv files */
template<class M, class V, class E, class P> inline
void print_cells_csv(AbstractMesh<M,E,V,P> &m, const std::unordered_map<int, int> &tmp_labels_map,
                     const std::string filename, Profiler Prof) {
    Prof.push("Print Region Cells");
    std::ofstream fp;
    fp.open(filename.c_str());
    assert(fp.is_open());
    fp << std::fixed;
    fp.precision(2);

    fp << "# cell id, centroid x, centroid y, centroid z, region, subregion, field" << std::endl;
    for (uint pid=0; pid<m.num_polys(); ++pid) {
        int tmp_label = m.poly_data(pid).label;
        auto it = tmp_labels_map.find(tmp_label);
        int old_label = (it == tmp_labels_map.end()) ? tmp_label : it->second;
        vec3d c = m.poly_centroid(pid);
        fp << pid << ", " << c.x() << ", " << c.y() << ", " << c.z() << ", "
           << old_label << ", " << tmp_label << ", " << m.poly_data(pid).quality << std::endl;
    }
    fp.close();
    Prof.pop();
}

/**********************************************************************/

void shp_poly_add_vert(const vec3d &v, uint &count, std::vector<double> &x, std::vector<double> &y) {
    x[count] = v.x();
    y[count] = v.y();
    ++count;
}

/* print the vertices of the regions on separate shapefiles */
void print_regions_shp(Polygonmesh<> &m, const Statistics &stats,
                       const std::unordered_map<int, std::vector<uint>> &labels_polys_map,
                       const std::string filename, Profiler Prof) {
    Prof.push("Print Region Shapefiles");
    // create the shapefile
    SHPHandle hSHP = SHPCreate((filename + ".shp").c_str(), SHPT_POLYGON);
    if (hSHP == nullptr) {
        std::cerr << "Error creating shapefile." << std::endl;
        return;
    }

    // create the attributes file
    DBFHandle hDBF = DBFCreate((filename + ".dbf").c_str());
    if (hDBF == nullptr) {
        std::cerr << "Error creating DBF file." << std::endl;
        SHPClose(hSHP);
        return;
    }
    DBFAddField(hDBF, "Label",      FTInteger, 50,  0);
    DBFAddField(hDBF, "N. Cells",   FTInteger, 1e9, 0);
    DBFAddField(hDBF, "Point",      FTString,  1e9, 0);
    DBFAddField(hDBF, "Size",       FTInteger, 1e9, 0);
    DBFAddField(hDBF, "F Mean",     FTDouble,  1e4, 3);
    DBFAddField(hDBF, "F StDev",    FTDouble,  1e4, 3);
    DBFAddField(hDBF, "F Min",      FTDouble,  1e4, 3);
    DBFAddField(hDBF, "F Max",      FTDouble,  1e4, 3);
    DBFAddField(hDBF, "Isoval Min", FTDouble,  1e4, 3);
    DBFAddField(hDBF, "Isoval Max", FTDouble,  1e4, 3);

    uint i = 0;
    for (auto &l : labels_polys_map) {
        std::vector<std::vector<uint>> verts_loops;
        extract_region(m, l.first, verts_loops);

        // find the outer loop (the connected component containing more vertices)
        // NOTE: it would be safer to pick the one with the greatest area
        auto it = std::max_element(verts_loops.begin(), verts_loops.end(),
                                   [](const std::vector<uint>& a, const std::vector<uint>& b) {
                                       return a.size() < b.size();
                                   });
        size_t index = std::distance(verts_loops.begin(), it);
        std::vector<uint> outer_loop = verts_loops.at(index);
        // remove outer_loop from the list
        if (it != verts_loops.end()) {
            verts_loops.erase(it);
        }

        // total number of vertices (outer_loop + inner_loops)
        int total_loops = 1 + verts_loops.size(); // +1 for outer_loop
        std::vector<int> loop_start_indices(total_loops);
        int total_vertices = 0;

        loop_start_indices.front() = total_vertices;
        total_vertices += outer_loop.size() + 1; // +1 for closing the loop

        for (int i=0; i<verts_loops.size(); ++i) {
            loop_start_indices.at(i+1) = total_vertices;
            total_vertices += verts_loops.at(i).size() + 1; // +1 for closing the loop
        }

        std::vector<double> x(total_vertices), y(total_vertices);
        // for windows users:
        // double *x = new double(n);
        // double *y = new double(n);
        uint count = 0;

        // outer_loop
        REVERSE_VEC(outer_loop);
        for (uint vid : outer_loop) {
            shp_poly_add_vert(m.vert(vid), count, x, y);
        }
        vec3d v = m.vert(outer_loop.front());
        shp_poly_add_vert(v, count, x, y);

        // inner_loop
        for (auto &loop : verts_loops) {
            REVERSE_VEC(loop);
            for (uint vid : loop) {
                shp_poly_add_vert(m.vert(vid), count, x, y);
            }
            vec3d v = m.vert(loop.front());
            shp_poly_add_vert(v, count, x, y);
        }

        // create the shape object
        SHPObject* obj = SHPCreateObject(SHPT_POLYGON, -1, total_loops, loop_start_indices.data(),
                                         nullptr, total_vertices, x.data(), y.data(), nullptr, nullptr);

        // write polygon on file with attributes
        int objId = SHPWriteObject(hSHP, -1, obj);

        if (!DBFWriteIntegerAttribute(hDBF, objId, 0, stats.label_vec[i])) {
            std::cerr << "Error writing Label attribute to DBF file." << std::endl;
        }
        if (!DBFWriteIntegerAttribute(hDBF, objId, 1, stats.card_vec[i])) {
            std::cerr << "Error writing N. Cells attribute to DBF file." << std::endl;
        }
        vec3d p = stats.pos_vec[i];
        std::string ss = std::to_string(p.x()) + ", " + std::to_string(p.y()) + ", " + std::to_string(p.z());
        if (!DBFWriteStringAttribute(hDBF, objId, 2, ss.c_str())) {
            std::cerr << "Error writing Point attribute to DBF file." << std::endl;
        }
        if (!DBFWriteIntegerAttribute(hDBF, objId, 3, stats.size_vec[i])) {
            std::cerr << "Error writing Size attribute to DBF file." << std::endl;
        }
        if (!DBFWriteDoubleAttribute(hDBF, objId, 4, stats.mean_vec[i])) {
            std::cerr << "Error writing F Mean attribute to DBF file." << std::endl;
        }
        if (!DBFWriteDoubleAttribute(hDBF, objId, 5, stats.stdev_vec[i])) {
            std::cerr << "Error writing F StDev attribute to DBF file." << std::endl;
        }
        if (!DBFWriteDoubleAttribute(hDBF, objId, 6, stats.minmax_vec[i].first)) {
            std::cerr << "Error writing F Min attribute to DBF file." << std::endl;
        }
        if (!DBFWriteDoubleAttribute(hDBF, objId, 7, stats.minmax_vec[i].second)) {
            std::cerr << "Error writing F Max attribute to DBF file." << std::endl;
        }
        if (!DBFWriteDoubleAttribute(hDBF, objId, 8, stats.iso_vec[i].first)) {
            std::cerr << "Error writing Isoval Min attribute to DBF file." << std::endl;
        }
        if (!DBFWriteDoubleAttribute(hDBF, objId, 9, stats.iso_vec[i].second)) {
            std::cerr << "Error writing Isoval Max attribute to DBF file." << std::endl;
        }
        SHPDestroyObject(obj);
        ++i;
    }
    SHPClose(hSHP);
    DBFClose(hDBF);
    Prof.pop();
}

/**********************************************************************/

/* print the regions as separate meshes, grouped for label */
void print_regions_off(Polygonmesh<> &m, const std::unordered_map<int, int> &tmp_labels_map,
                       const int n_labels, const std::string base_path, Profiler Prof) {
    Prof.push("Print Region Meshes");
    std::string dir = base_path + "/regions";
    open_directory(dir);

    for (uint l=0; l<n_labels; ++l) {
        // extract the mesh containing only the cells with label *l*
        Polygonmesh<> sub_m;
        std::unordered_map<uint,uint> m2subm_vmap, subm2m_vmap;
        export_cluster(m, l, sub_m, m2subm_vmap, subm2m_vmap);
        if (sub_m.num_polys() == 0)
            continue;
        std::string path = dir + "/region_" + std::to_string(l) + ".off";
        sub_m.save(path.c_str());

        // recover verts and cells data from *m* and print them in csv files
        for (uint pid_new=0; pid_new<sub_m.num_polys(); ++pid_new) {
            std::vector<uint> vids_new = sub_m.poly_verts_id(pid_new);
            std::vector<uint> vids_old;
            for (uint vid : vids_new) {
                vids_old.push_back(subm2m_vmap[vid]);
            }
            int pid_old = m.poly_id(vids_old);
            assert(pid_old>=0 && "print_regions_off: pid_old not found!");
            sub_m.poly_data(pid_new) = m.poly_data(pid_old);
        }
        std::string path2 = dir + "/region_" + std::to_string(l);
        print_verts_csv(sub_m, path2 + "_verts_data.csv", Prof);
        print_cells_csv(sub_m, tmp_labels_map, path2 + "_cells_data.csv", Prof);
    }
    Prof.pop();
}

/* print the regions as separate meshes, grouped for label */
void print_regions_off(Polyhedralmesh<> &m, const std::unordered_map<int, int> &tmp_labels_map,
                       const int n_labels, const std::string base_path, Profiler Prof) {
    Prof.push("Print Region Meshes");
    std::string dir = base_path + "/regions";
    open_directory(dir);

    for (uint l=0; l<n_labels; ++l) {
        // extract the mesh containing only the cells with label *l*
        Polyhedralmesh<> sub_m;
        std::unordered_map<uint,uint> m2subm_vmap, subm2m_vmap;
        export_cluster(m, l, sub_m, m2subm_vmap, subm2m_vmap);
        if (sub_m.num_polys() == 0)
            continue;
        std::string path = dir + "/region_" + std::to_string(l) + ".off";
        sub_m.save(path.c_str());

        // recover verts and cells data from *m* and print them in csv files
        for (uint pid_new=0; pid_new<sub_m.num_polys(); ++pid_new) {
            std::vector<uint> vids_new = sub_m.poly_verts_id(pid_new);
            std::vector<uint> vids_old;
            for (uint vid : vids_new) {
                vids_old.push_back(subm2m_vmap[vid]);
            }
            int pid_old = m.poly_id(vids_old);
            assert(pid_old>=0 && "print_regions_off: pid_old not found!");
            sub_m.poly_data(pid_new) = m.poly_data(pid_old);
        }
        std::string path2 = dir + "/region_" + std::to_string(l);
        print_verts_csv(sub_m, path2 + "_verts_data.csv", Prof);
        print_cells_csv(sub_m, tmp_labels_map, path2 + "_cells_data.csv", Prof);
    }
    Prof.pop();
}

/**********************************************************************/

/* print the regions boundaries as separate meshes, grouped for label */
void print_regions_bnd(Polygonmesh<> &m, const std::unordered_map<int, std::vector<uint>> &labels_polys_map,
                       const std::vector<int> &label_vec, const std::string base_path, Profiler Prof) {
    Prof.push("Print Region Meshes");
    std::string dir = base_path + "/regions_bnd";
    open_directory(dir);
    //
    for (uint label : label_vec) {
        std::vector<std::vector<uint>> verts;
        extract_region(m, label, verts);

        Polygonmesh<> sub_m(m.vector_verts());
        for (std::vector<uint> poly : verts) {
            sub_m.poly_add(poly);
        }

        std::string path = dir + "/region_" + std::to_string(label) + ".off";
        sub_m.save(path.c_str());
    }
    Prof.pop();
}

/**********************************************************************/

/* print the cell labels in a csv file */
template<class M, class V, class E, class P> inline
void print_labels(AbstractMesh<M,E,V,P> &m, const std::string output_path)
{
    std::string labels_path = output_path + "/labels.csv";
    std::ofstream fp;
    fp.open(labels_path.c_str());
    assert(fp.is_open());
    fp << "LABEL_FIELD " << m.num_polys() << "\n";
    for (uint pid = 0; pid < m.num_polys(); ++pid) {
        fp << m.poly_data(pid).label << std::endl;
    }
    fp.close();
}

/**********************************************************************/

/* save the final mesh */
template<class M, class V, class E, class P> inline
void save_mesh(AbstractMesh<M,E,V,P> &m,
               const std::string output_path, const std::string MESH_FORMAT)
{
    if (MESH_FORMAT == "-") {
        MESH_FORMAT = m.mesh_is_surface() ? "obj" : "vtk";
    }
    std::string mesh_path = output_path + "/mesh_out." + MESH_FORMAT;
    m.save(mesh_path.c_str());
}

#endif // AUXILIARY_H
