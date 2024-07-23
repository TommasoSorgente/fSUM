#ifndef STATISTICS_H
#define STATISTICS_H

/* class for printing geometrical and statistical information about the regions of the partition */
 
#include <cinolib/export_cluster.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/profiler.h>
#include <cinolib/connected_components.h>

#include <iostream>
#include <filesystem>
#include <shapefil.h>
#include <cstddef>

namespace fs = std::filesystem;
using namespace cinolib;

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

class Statistics {

private:
  std::string base_path;
  std::vector<double> perc_vals;
  std::vector<double> isovals;
  std::unordered_map<int, std::vector<uint>> labels_polys_map;
  Profiler Prof;

public:
  Statistics(std::string filename, std::vector<double> perc, std::vector<double> vals,
             std::unordered_map<int, std::vector<uint>> map) {
    size_t last = filename.find_last_of('/');
    std::string name_tmp = filename.substr(0, last);
    size_t first = name_tmp.find_last_of('/');
    std::string name = name_tmp.substr(first+1);
    base_path = std::string(HOME_PATH) + "out/" + name;
    open_directory(base_path);

    perc_vals = perc;
    isovals = vals;
    labels_polys_map = map;
  }
  ~Statistics() {}

  template <class M, class V, class E, class P>
  inline void compute_stats(AbstractMesh<M, E, V, P> &m);

  template <class M, class V, class E, class P>
  inline void print_stats(AbstractMesh<M, E, V, P> &m);

  void extract_region(Polygonmesh<> &m, const uint label,
                      std::vector<std::vector<vec3d>> &verts);
  void print_regions_csv(Polygonmesh<> &m);
  void print_regions_shp(Polygonmesh<> &m);

  // vectors for storing information about each region:
  std::vector<int> father_vec;   // father of the region
  std::vector<int> label_vec;    // label of the region
  std::vector<uint> card_vec;    // number of cells in the region
  std::vector<vec3d> pos_vec;    // random point in the region
  std::vector<double> size_vec;  // area of the region
  std::vector<double> mean_vec;  // mean field value in the region
  std::vector<double> stdev_vec; // standard deviation of the field
  std::vector<std::pair<double, double>> minmax_vec; // field extremities in the region
  std::vector<std::pair<double, double>> iso_vec;    // isovalues of the region
};

/**********************************************************************/

template <class M, class V, class E, class P>
inline void Statistics::compute_stats(AbstractMesh<M, E, V, P> &m) {
  Prof.push("Compute Statistics");
  for (auto &l : labels_polys_map) {
    int n = l.second.size();
    uint cid = l.second.front();
    //
    father_vec.push_back(l.first); // TBD
    label_vec.push_back(l.first);
    card_vec.push_back(n);
    pos_vec.push_back(m.poly_centroid(cid));
    //
    double mean = 0., sqmean = 0., stdev = 0., size = 0.;
    double fmin = DBL_MAX, fmax = DBL_MIN;
    for (uint pid : l.second) {
      double val = m.poly_data(pid).quality;
      mean += val;
      sqmean += val * val;
      size += m.poly_mass(pid);
      fmin = std::min(fmin, val);
      fmax = std::max(fmax, val);
    }
    mean /= n;
    stdev = sqrt(sqmean / n - mean * mean);
    //
    mean_vec.push_back(mean);
    stdev_vec.push_back(stdev);
    size_vec.push_back(size);
    minmax_vec.push_back(std::pair<double, double>(fmin, fmax));

    double val = m.poly_data(cid).quality;
    for (uint lid = 0; lid < isovals.size() - 1; ++lid) {
      if (isovals.at(lid) <= val && val <= isovals.at(lid + 1)) {
        iso_vec.push_back(std::pair<double, double>(perc_vals.at(lid),
                                                    perc_vals.at(lid + 1)));
        break;
      }
    }
  }
  Prof.pop();
}

/**********************************************************************/

template <class M, class V, class E, class P>
inline void Statistics::print_stats(AbstractMesh<M, E, V, P> &m) {
  Prof.push("Print Statistics");
  std::string path = base_path + "/stats";
  open_directory(path);
  //
  for (uint i = 0; i < labels_polys_map.size(); ++i) {
    std::string filename = path + "/region_" + std::to_string(i) + ".csv";
    std::ofstream fp;
    fp.open(filename.c_str());
    assert(fp.is_open());
    //
    fp << std::fixed;
    fp.precision(0);
    fp << "FATHER: " << father_vec[i] << std::endl;
    fp << "LABEL: " << label_vec[i] << std::endl;
    fp << "ISO_REG: " << iso_vec[i].first << " - " << iso_vec[i].second << std::endl;
    fp << "#CELLS: " << card_vec[i] << std::endl;
    fp << "POINT: " << pos_vec[i].x() << " " << pos_vec[i].y() << " " << pos_vec[i].z() << std::endl;
    fp.precision(2);
    fp << "SIZE: " << size_vec[i] << std::endl;
    fp << "MEAN: " << mean_vec[i] << std::endl;
    fp << "STDEV: " << stdev_vec[i] << std::endl;
    fp << "MIN: " << minmax_vec[i].first << std::endl;
    fp << "MAX: " << minmax_vec[i].second << std::endl;

    // fp << "\nFIELD VALUES:\n";
    // std::vector<uint> polys = labels_polys_map.at(label_vec[i]);
    // for (uint pid : polys) {
    //     fp << m.poly_data(pid).quality << std::endl;
    // }
    fp.close();
  }
  Prof.pop();
}

/**********************************************************************/

void Statistics::extract_region(Polygonmesh<> &m, const uint label,
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

/**********************************************************************/

void Statistics::print_regions_csv(Polygonmesh<> &m) {
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

void Statistics::print_regions_shp(Polygonmesh<> &m) {
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

#endif // STATISTICS_H
