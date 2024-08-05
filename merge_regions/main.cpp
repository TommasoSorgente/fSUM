#include <iostream>
#include <filesystem>

#include <cinolib/meshes/meshes.h>
#include <cinolib/profiler.h>
#include <cinolib/merge_meshes_at_coincident_vertices.h>

#include "../src/auxiliary.h"
#include "../src/statistics.h"

using namespace cinolib;
namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
    Profiler    prof;

    std::vector<std::string> regions{"cr_F5TERRE",
                                     "cr_FMAGRA",
                                     "cr_FPETRONIO",
                                     "cr_FPORTOFINO",
                                     "cr_FENTELLA",
                                     "cr_FAVETO",
                                     "cr_FPADANO",
                                     "cr_FBISAGNO",
                                     "cr_FPOLCEVERA",
                                     "cr_FARENZANO",
                                     "cr_FSASSELLO",
                                     "cr_FBORMIDE",
                                     "cr_FSAVONESE",
                                     "cr_FIMPERIESE"};

    std::string field_path = std::string(HOME_PATH) + "../data/Liguria_grid/cr_mean_global.csv";
    std::vector<double> field_global;
    std::ifstream fp(field_path.c_str());
    assert(fp.is_open());
    double f;
    while (fp >> f) {
        field_global.push_back(f);
    }
    fp.close();

    uint n_labels = 4;
    std::vector<double> perc_vals = {0, 15, 45, 90, 100};
    std::vector<double> isovals;
    for (int val : perc_vals) {
        isovals.push_back(percentile(field_global, val));
    }

    for (uint label=0; label<n_labels; ++label) {
        prof.push("Processing label " + std::to_string(label));
        uint count = 0;
        Polygonmesh<> m;

        for (uint rid=0; rid<regions.size(); ++rid) {
            std::string region_mesh  = std::string(HOME_PATH) + "../out/" +
                                       regions.at(rid) + "/regions_off/region_" + std::to_string(label) + ".off";
            std::string region_field = std::string(HOME_PATH) + "../data/Liguria_grid/" +
                                       regions.at(rid) + "/cr_mean.csv";
            if (!fs::exists(region_mesh)) {
                std::cout << "Warning: region " << region_mesh << " does not exist!" << std::endl;
                continue;
            }

            Polygonmesh<> m2(region_mesh.c_str());

            // load the field from file
            std::map<uint,double> field;
            std::ifstream fp(region_field.c_str());
            assert(fp.is_open());
            double f;
            int _pid = 0;
            while (fp >> f) {
                field[_pid] = f;
                ++_pid;
            }
            fp.close();

            // assert(field.size() == m2.num_polys());
            for (uint pid=0; pid<m2.num_polys(); ++pid) {
                m2.poly_data(pid).quality = field.at(pid);
                m2.poly_data(pid).label = label;
            }

            Polygonmesh<> res;

            merge_meshes_at_coincident_vertices(m, m2, res);
            for (uint pid=0; pid<m.num_polys(); ++pid) {
                res.poly_data(pid) = m.poly_data(pid);
            }
            int offset = m.num_polys();
            for (uint pid=0; pid<m2.num_polys(); ++pid) {
                res.poly_data(offset + pid) = m2.poly_data(pid);
            }

            std::unordered_map<int, std::vector<uint>> labels_polys_map;
            update_labels_polys_map(res, labels_polys_map);
            std::vector<int> label_vec;
            for (auto &l : labels_polys_map) {
                label_vec.push_back(l.first);
            }
            std::string output_path = std::string(HOME_PATH) + "out/label_" + std::to_string(label);

            Statistics stats(output_path, perc_vals, isovals, labels_polys_map, prof);
            stats.compute_stats(res);
            stats.print_stats(res);

            print_verts_csv(res, labels_polys_map, label_vec, output_path, prof);
            print_cells_csv(res, labels_polys_map, label_vec, output_path, prof);
            print_regions_shp(res, labels_polys_map, label_vec, output_path, prof);
            print_regions_off(res, labels_polys_map, label_vec, output_path, prof);
            // print_regions_bnd(res, labels_polys_map, label_vec, output_path, prof);

            print_labels(res, output_path);
            // save_mesh(res, output_path, Par.get_OUT_FORMAT());

            m = res;
            ++count;
        }
        prof.pop(true, " regions found: " + std::to_string(count));
        std::cout << std::endl;
    }

    return 0;
}
