#ifndef AUXILIARY_H
#define AUXILIARY_H
 
#include <filesystem>

#include <cinolib/meshes/meshes.h>
#include <cinolib/export_cluster.h>
#include <cinolib/profiler.h>

#include "statistics.h"

using namespace cinolib;
namespace fs = std::filesystem;

/* Auxiliary functions */

/**********************************************************************/

// compute the isovalue corresponding to the percentile in *data*
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
    double isovalue = sortedData[index] + interpolation * (sortedData[index+1] - sortedData[index]);
    return isovalue;
}

// compute the percentile corresponding to the isovalue in *data*
double inv_percentile(const std::vector<double> &data, const double isovalue)
{
    assert(!data.empty());
    std::vector<double> sortedData = data;
    sort(sortedData.begin(), sortedData.end());
    int n = sortedData.size()-1;

    auto it = std::lower_bound(sortedData.begin(), sortedData.end(), isovalue);
    int index = it - sortedData.begin();
    double percentile = index * 100 / n;
    assert(0. <= percentile && percentile <= 100.);
    return percentile;
}

/**********************************************************************/

// compute the mean value of *field* in the 1-ring of poly *pid*
template<class M, class E, class V, class P> inline
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
template<class M, class E, class V, class P> inline
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
template<class M, class E, class V, class P> inline
void mark_edges(AbstractMesh<M,E,V,P> &m)
{
    for (uint eid=0; eid<m.num_edges(); ++eid) {
        m.edge_data(eid).flags[MARKED] = false;
        std::vector<uint> polys = m.adj_e2p(eid);
        if (m.poly_data(polys.front()).label != m.poly_data(polys.back()).label) {
            m.edge_data(eid).flags[MARKED] = true;
        }
    }
}

void mark_faces(Polyhedralmesh<> &m)
{
    for (uint fid=0; fid<m.num_faces(); ++fid) {
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
    }
}

/**********************************************************************/

/* create a map <label, all polys with that label> */
template<class M, class E, class V, class P> inline
void update_regions_map(AbstractMesh<M,E,V,P> &m, std::unordered_map<int, std::vector<uint>> &map)
{
    map.clear();
    for (uint pid = 0; pid < m.num_polys(); ++pid) {
        map[m.poly_data(pid).label].push_back(pid);
    }
}

/* apply the labels in *map* to the mesh elements */
template<class M, class E, class V, class P> inline
void apply_labels_to_mesh(AbstractMesh<M,E,V,P> &m, std::unordered_map<int, std::vector<uint>> &map)
{
    for (auto l : map) {
        for (uint pid : l.second) {
            m.poly_data(pid).label = l.first;
        }
    }
}

/**********************************************************************/

template<class M, class E, class V, class P> inline
double poly_misclassification(AbstractMesh<M,E,V,P> &m, const uint pid, const std::pair<double,double> &isovals, const std::pair<double,double> &sigma)
{
    double f = m.poly_data(pid).quality;
    if (sigma.second < isovals.first) {
        return isovals.first - f;
    }
    else if (sigma.first > isovals.second) {
        return f - isovals.second;
    }
    return 0.;
}

/**********************************************************************/

double mesh_shape_coefficient(Polygonmesh<> &m, Polygonmesh<> &sub_m) {
    double perimeter = 0.;
    for (std::pair<uint,uint> e : sub_m.get_boundary_edges()) {
        perimeter += sub_m.vert(e.first).dist(sub_m.vert(e.second));
    }
    AABB bbox = sub_m.bbox();
    double bbox_perimeter = 2. * (bbox.delta_x() + bbox.delta_y());
    double area =sub_m.mesh_area();
    double entr = area / m.mesh_area();

    double compactness = sqrt(area) / perimeter;
    double smoothness  = bbox_perimeter / perimeter;
    double entropy     = -entr * log(entr);
    return (compactness + smoothness + entropy) / 3.;
}

double mesh_shape_coefficient(Polyhedralmesh<> &m, Polyhedralmesh<> &sub_m) {
    double srf_area = sub_m.mesh_srf_area();
    AABB bbox = sub_m.bbox();
    double bbox_area = 2. * (bbox.delta_x() * bbox.delta_y() +
                             bbox.delta_x() * bbox.delta_z() +
                             bbox.delta_y() * bbox.delta_z());
    double volume = sub_m.mesh_volume();
    double entr   = volume / m.mesh_volume();

    double compactness = cbrt(volume) / srf_area;
    double smoothness  = bbox_area / srf_area;
    double entropy     = entr * log(entr);
    return (compactness + smoothness + entropy) / 3.;
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

/* print scalar field */
template<class M, class E, class V, class P> inline
void print_field_csv(AbstractMesh<M,E,V,P> &m, const std::string filename, Profiler Prof) {
    Prof.push("Print Field");
    std::ofstream fp(filename.c_str());
    assert(fp.is_open());
    fp << std::fixed;

    fp << "# x, y, z, field" << std::endl;
    for (uint pid=0; pid<m.num_polys(); ++pid) {
        vec3d c = m.poly_centroid(pid);
        fp.precision(2);
        fp << c.x() << ", " << c.y() << ", " << c.z() << ", ";
        fp.precision(6);
        fp << m.poly_data(pid).quality << std::endl;
    }
    fp.close();
    Prof.pop();
}

/**********************************************************************/

/* print misclassifications */
template<class M, class E, class V, class P> inline
void print_misclassifications(AbstractMesh<M,E,V,P> &m, const std::vector<double> &isovals,
                              const std::unordered_map<int, int> &subregion_region_map,
                              const std::string filename, Profiler Prof) {
    Prof.push("Print Misclassifications");
    std::ofstream fp(filename.c_str());
    assert(fp.is_open());

    fp << "# pid, centroid x, centroid y, centroid z, field, isovalue_min, isovalue_max, misclassification" << std::endl;
    fp << std::fixed;

    for(uint pid=0; pid<m.num_polys(); ++pid) {
        int tmp_label = m.poly_data(pid).label;
        auto it = subregion_region_map.find(tmp_label);
        int l = (it == subregion_region_map.end()) ? tmp_label : it->second;
        double iso_min = isovals.at(l);
        double iso_max = isovals.at(l+1);

        double misclass = poly_misclassification(m, pid, {iso_min,iso_max}, {0.,0.});
        vec3d c = m.poly_centroid(pid);
        if (misclass > 0) {
            fp.precision(2);
            fp << pid << ", " << c.x() << ", " << c.y() << ", " << c.z() << ", ";
            fp.precision(6);
            fp << m.poly_data(pid).quality << ", " << iso_min << ", " << iso_max << ", " << misclass << std::endl;
        }
    }
    fp.close();
    Prof.pop();
}

/**********************************************************************/

/* print isovalues info */
template<class M, class E, class V, class P> inline
void print_isovals_csv(AbstractMesh<M,E,V,P> &m, const std::vector<double> &percentiles, const std::vector<double> &isovals, const std::string filename, Profiler Prof) {
    Prof.push("Print Isovalues");
    std::ofstream fp(filename.c_str());
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
template<class M, class E, class V, class P> inline
void print_verts_csv(AbstractMesh<M,E,V,P> &m, const std::string filename, Profiler Prof) {
    Prof.push("Print Region Vertices");
    std::ofstream fp(filename.c_str());
    assert(fp.is_open());
    fp << std::fixed;

    fp << "# vertex id, x, y, z, field" << std::endl;
    for (uint vid=0; vid<m.num_verts(); ++vid) {
        vec3d v = m.vert(vid);
        fp.precision(2);
        fp << vid << ", " << v.x() << ", " << v.y() << ", " << v.z() << ", ";
        fp.precision(6);
        fp << m.vert_data(vid).quality << std::endl;
    }
    fp.close();
    Prof.pop();
}

/**********************************************************************/

/* print the cells of the regions on separate csv files */
template<class M, class E, class V, class P> inline
void print_cells_csv(AbstractMesh<M,E,V,P> &m, const std::unordered_map<int, int> &subregion_region_map,
                     const std::string filename, Profiler Prof) {
    Prof.push("Print Region Cells");
    std::ofstream fp(filename.c_str());
    assert(fp.is_open());
    fp << std::fixed;

    fp << "# cell id, centroid x, centroid y, centroid z, region, subregion, field" << std::endl;
    for (uint pid=0; pid<m.num_polys(); ++pid) {
        int tmp_label = m.poly_data(pid).label;
        auto it = subregion_region_map.find(tmp_label);
        int old_label = (it == subregion_region_map.end()) ? tmp_label : it->second;
        vec3d c = m.poly_centroid(pid);
        fp.precision(2);
        fp << pid << ", " << c.x() << ", " << c.y() << ", " << c.z()
           << ", " << old_label << ", " << tmp_label << ", ";
        fp.precision(6);
        fp << m.poly_data(pid).quality << std::endl;
    }
    fp.close();
    Prof.pop();
}

/**********************************************************************/

/* print the regions as separate meshes, grouped for label */
void print_regions(Polygonmesh<> &m, const std::unordered_map<int, int> &subregion_region_map, std::unordered_map<int, std::vector<uint>> polys_in_region,
Statistics &stats, const std::string base_path, Profiler Prof) {
    Prof.push("Print Region Meshes");
    open_directory(base_path);

    for (auto l : polys_in_region) {        
        // extract the mesh containing only the cells with label *l*
        Polygonmesh<> sub_m;
        std::unordered_map<uint,uint> m2subm_vmap, subm2m_vmap;
        export_cluster(m, l.first, sub_m, m2subm_vmap, subm2m_vmap);
        if (sub_m.num_polys() == 0)
            continue;
        std::string path = base_path + "/region_" + std::to_string(l.first) + ".obj";
        sub_m.save(path.c_str());

        // recover verts and cells data from *m* and print them in csv files
        for (uint pid_new=0; pid_new<sub_m.num_polys(); ++pid_new) {
            std::vector<uint> vids_new = sub_m.poly_verts_id(pid_new);
            std::vector<uint> vids_old;
            for (uint vid : vids_new) {
                vids_old.push_back(subm2m_vmap[vid]);
            }
            int pid_old = m.poly_id(vids_old);
            assert(pid_old>=0 && "print_regions: pid_old not found!");
            sub_m.poly_data(pid_new) = m.poly_data(pid_old);
        }
        std::string path2 = base_path + "/region_" + std::to_string(l.first);
        print_verts_csv(sub_m, path2 + "_verts_data.csv", Prof);
        print_cells_csv(sub_m, subregion_region_map, path2 + "_cells_data.csv", Prof);

        std::unordered_map<int, std::vector<uint>> polys_in_region2;
        update_regions_map(sub_m, polys_in_region2);
        stats.clear_stats();
        stats.compute_stats(sub_m, polys_in_region2);
        stats.print_stats(sub_m, base_path);
    }
    Prof.pop();
}

/* print the regions as separate meshes, grouped for label */
void print_regions(Polyhedralmesh<> &m, const std::unordered_map<int, int> &subregion_region_map, std::unordered_map<int, std::vector<uint>> polys_in_region,
                   Statistics &stats, const std::string base_path, Profiler Prof) {
    Prof.push("Print Region Meshes");
    open_directory(base_path);

    for (auto l : polys_in_region) {
        // extract the mesh containing only the cells with label *l*
        Polyhedralmesh<> sub_m;
        std::unordered_map<uint,uint> m2subm_vmap, subm2m_vmap;
        export_cluster(m, l.first, sub_m, m2subm_vmap, subm2m_vmap);
        if (sub_m.num_polys() == 0)
            continue;
        std::string path = base_path + "/region_" + std::to_string(l.first) + ".vtk";
        sub_m.save(path.c_str());

        // recover verts and cells data from *m* and print them in csv files
        for (uint pid_new=0; pid_new<sub_m.num_polys(); ++pid_new) {
            std::vector<uint> vids_new = sub_m.poly_verts_id(pid_new);
            std::vector<uint> vids_old;
            for (uint vid : vids_new) {
                vids_old.push_back(subm2m_vmap[vid]);
            }
            int pid_old = m.poly_id(vids_old);
            assert(pid_old>=0 && "print_regions: pid_old not found!");
            sub_m.poly_data(pid_new) = m.poly_data(pid_old);
        }
        std::string path2 = base_path + "/region_" + std::to_string(l.first);
        print_verts_csv(sub_m, path2 + "_verts_data.csv", Prof);
        print_cells_csv(sub_m, subregion_region_map, path2 + "_cells_data.csv", Prof);

        std::unordered_map<int, std::vector<uint>> polys_in_region2;
        update_regions_map(sub_m, polys_in_region2);
        stats.clear_stats();
        stats.compute_stats(sub_m, polys_in_region2);
        stats.print_stats(sub_m, base_path);
    }
    Prof.pop();
}

/**********************************************************************/

// /* print the regions boundaries as separate meshes, grouped for label */
// void print_regions_bnd(Polygonmesh<> &m, const std::unordered_map<int, std::vector<uint>> &polys_in_region,
//                        const std::vector<int> &label_vec, const std::string base_path, Profiler Prof) {
//     Prof.push("Print Region Meshes");
//     std::string dir = base_path + "/regions_bnd";
//     open_directory(dir);
//     //
//     for (uint label : label_vec) {
//         std::vector<std::vector<uint>> verts;
//         extract_region(m, label, verts);

//         Polygonmesh<> sub_m(m.vector_verts());
//         for (std::vector<uint> poly : verts) {
//             sub_m.poly_add(poly);
//         }

//         std::string path = dir + "/region_" + std::to_string(label) + ".off";
//         sub_m.save(path.c_str());
//     }
//     Prof.pop();
// }

/**********************************************************************/

/* print the cell labels in a csv file */
template<class M, class E, class V, class P> inline
void print_labels(AbstractMesh<M,E,V,P> &m, const std::string output_path)
{
    std::string labels_path = output_path + "/labels.csv";
    std::ofstream fp(labels_path.c_str());
    assert(fp.is_open());
    fp << "LABEL_FIELD " << m.num_polys() << "\n";
    for (uint pid = 0; pid < m.num_polys(); ++pid) {
        fp << m.poly_data(pid).label << std::endl;
    }
    fp.close();
}

/**********************************************************************/

/* save the final mesh */
template<class M, class E, class V, class P> inline
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
