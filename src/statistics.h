#ifndef STATISTICS_H
#define STATISTICS_H

#include <cinolib/export_cluster.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/profiler.h>
#include <cinolib/connected_components.h>

#include "../src/auxiliary.h"

using namespace cinolib;

/* class for printing geometrical and statistical information about the regions of the partition */
class Statistics {

private:
  std::string base_path;
  std::vector<double> perc_vals;
  std::vector<double> isovals;
  std::unordered_map<int, std::vector<uint>> labels_polys_map;
  Profiler Prof;

public:
  Statistics(std::string output_path, std::vector<double> perc, std::vector<double> vals,
             std::unordered_map<int, std::vector<uint>> map, Profiler p) {
    base_path = output_path;
    perc_vals = perc;
    isovals = vals;
    labels_polys_map = map;
    Prof = p;
  }
  ~Statistics() {}

  template <class M, class V, class E, class P>
  inline void compute_stats(AbstractMesh<M, E, V, P> &m);

  template <class M, class V, class E, class P>
  inline void print_stats(AbstractMesh<M, E, V, P> &m);

  // vectors for storing information about each region:
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

#endif // STATISTICS_H
