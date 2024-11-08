#ifndef STATISTICS_H
#define STATISTICS_H

#include <cinolib/export_cluster.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/profiler.h>
#include <cinolib/connected_components.h>

using namespace cinolib;

/* class for printing geometrical and statistical information about the regions of the partition */
class Statistics {

private:
  std::vector<double> percvals, isovals;
  std::unordered_map<int, std::vector<uint>> polys_in_region;
  Profiler Prof;

  void load_sigma(const std::string path, std::vector<double> &field) {
      if (path == "") return;
      std::ifstream fp(path.c_str());
      assert(fp.is_open());
      double f;
      while (fp >> f) {
          field.push_back(f);
      }
      fp.close();
  }

public:
  Statistics() {}
  ~Statistics() {}

  void init(const std::vector<double> &perc, const std::vector<double> &vals, Profiler p,
            const std::string sigma_m="", const std::string sigma_p="") {
      percvals = perc;
      isovals = vals;
      Prof = p;
      polys_in_region.clear();
      clear_stats();
  }

  void clear_stats();

  template<class M, class V, class E, class P> inline
  void compute_stats(AbstractMesh<M,V,E,P> &m, std::unordered_map<int, std::vector<uint>> map);

  template<class M, class V, class E, class P> inline
  void print_stats(AbstractMesh<M,V,E,P> &m, const std::string path);

  // vectors for storing information about each region:
  std::vector<int> label_vec;    // label of the region
  std::vector<uint> card_vec;    // number of cells in the region
  std::vector<vec3d> pos_vec;    // random point in the region
  std::vector<double> size_vec;  // area of the region
  std::vector<double> mean_vec;  // mean field value in the region
  std::vector<double> stdev_vec; // standard deviation of the field
  std::vector<std::pair<double, double>> minmax_vec; // field extremities in the region
  std::vector<std::pair<double, double>> perc_vec;   // percentiles of the region
  std::vector<std::pair<double, double>> iso_vec;    // isovalues of the region
};

/**********************************************************************/

#include "auxiliary.h"

void Statistics::clear_stats() {
    label_vec.clear();
    card_vec.clear();
    pos_vec.clear();
    size_vec.clear();
    mean_vec.clear();
    stdev_vec.clear();
    minmax_vec.clear();
    perc_vec.clear();
    iso_vec.clear();
}

/**********************************************************************/

template<class M, class V, class E, class P> inline
void Statistics::compute_stats(AbstractMesh<M,V,E,P> &m, std::unordered_map<int, std::vector<uint>> map) {
    Prof.push("Compute Statistics");

    polys_in_region = map;
    for (auto &l : polys_in_region) {
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
            double val = m.poly_data(pid).fvalue;
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

        double val = mean; //m.poly_data(cid).fvalue;
        for (uint lid = 0; lid < isovals.size() - 1; ++lid) {
            if (isovals.at(lid) <= val && val <= isovals.at(lid + 1)) {
                perc_vec.push_back(std::pair<double, double>(percvals.at(lid),
                                                             percvals.at(lid + 1)));
                iso_vec.push_back(std::pair<double, double>(isovals.at(lid),
                                                            isovals.at(lid + 1)));
                break;
            }
        }
    }
    Prof.pop();
}

/**********************************************************************/

template<class M, class V, class E, class P> inline
void Statistics::print_stats(AbstractMesh<M,V,E,P> &m, const std::string path) {
  Prof.push("Print Statistics");

  for (uint i = 0; i < polys_in_region.size(); ++i) {
    std::string filename = path + "/region_" + std::to_string(label_vec[i]) + "_stats.csv";
    std::ofstream fp(filename.c_str());
    assert(fp.is_open());
    //
    fp << std::fixed;
    fp.precision(0);
    fp << "LABEL: " << label_vec[i] << std::endl;
    fp << "PERCENTILES: " << perc_vec[i].first << " - " << perc_vec[i].second << std::endl;
    fp.precision(4);
    fp << "ISOVALUES: " << iso_vec[i].first << " - " << iso_vec[i].second << std::endl;
    fp << "#CELLS: " << card_vec[i] << std::endl;
    fp.precision(2);
    fp << "POINT: " << pos_vec[i].x() << " " << pos_vec[i].y() << " " << pos_vec[i].z() << std::endl;
    fp << "SIZE: " << size_vec[i] << std::endl;
    fp << "MEAN: " << mean_vec[i] << std::endl;
    fp << "STDEV: " << stdev_vec[i] << std::endl;
    fp << "MIN: " << minmax_vec[i].first << std::endl;
    fp << "MAX: " << minmax_vec[i].second << std::endl;
    fp.close();
  }
  Prof.pop();
}

#endif // STATISTICS_H
