#ifndef SHAPEFILE_UTILS_H
#define SHAPEFILE_UTILS_H
 
#include <shapefil.h>

#include <cinolib/meshes/meshes.h>
#include <cinolib/export_cluster.h>
#include <cinolib/profiler.h>
#include <cinolib/connected_components.h>

#include "statistics.h"

using namespace cinolib;

/* Shapelib auxiliary functions */

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

void shp_poly_add_vert(const vec3d &v, uint &count, std::vector<double> &x, std::vector<double> &y) {
    x[count] = v.x();
    y[count] = v.y();
    ++count;
}

/**********************************************************************/

/* print the vertices of the regions on separate shapefiles */
void print_regions_shp(Polygonmesh<> &m, const Statistics &stats,
                       const std::unordered_map<int, std::vector<uint>> &polys_in_region,
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
    for (auto &l : polys_in_region) {
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

#endif // SHAPEFILE_UTILS_H
