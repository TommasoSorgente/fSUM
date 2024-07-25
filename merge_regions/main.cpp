#include <iostream>
#include <filesystem>

#include <cinolib/meshes/meshes.h>
#include <cinolib/profiler.h>
#include <cinolib/merge_meshes_at_coincident_vertices.h>

#include "../src/auxiliary.h"

using namespace cinolib;
namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
    std::string base_path = std::string(HOME_PATH) + "../out/";
    uint        n_labels  = 5;
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


    for (uint label = 0; label<n_labels; ++label) {

        prof.push("Processing label " + std::to_string(label));
        uint count = 0;
        Polygonmesh<> m;

        for (uint rid=0; rid<regions.size(); ++rid) {
            std::string region_name = base_path + regions.at(rid) + "/regions_off/region_"
                                      + std::to_string(label) + ".off";
            if (!fs::exists(region_name))
                continue;

            Polygonmesh<> m2(region_name.c_str());
            Polygonmesh<> res;

            merge_meshes_at_coincident_vertices(m, m2, res);
            for (uint pid=0; pid<m.num_polys(); ++pid) {
                res.poly_data(pid) = m.poly_data(pid);
            }
            for (uint pid=0; pid<m2.num_polys(); ++pid) {
                int offset = m.num_polys();
                res.poly_data(offset + pid) = m2.poly_data(pid);
            }

            std::unordered_map<int, std::vector<uint>> labels_polys_map;
            update_labels_polys_map(res, labels_polys_map);
            std::vector<int> label_vec;
            for (auto &l : labels_polys_map) {
                label_vec.push_back(l.first);
            }
            std::string output_path = std::string(HOME_PATH) + "/out/label_" + std::to_string(label);
            print_regions_shp(res, labels_polys_map, label_vec, output_path, prof);

            m = res;
            ++count;
        }
        prof.pop(true, "regions found: " + std::to_string(count));
        std::cout << std::endl;
    }

    return 0;
}
