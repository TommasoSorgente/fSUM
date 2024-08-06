#include <iostream>
#include <filesystem>

#include <cinolib/meshes/meshes.h>
#include <cinolib/profiler.h>
#include <cinolib/merge_meshes_at_coincident_vertices.h>

#include "../src/statistics.h"
#include "../src/auxiliary.h"

using namespace cinolib;
namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
    Profiler prof;

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
        prof.push("Processed label " + std::to_string(label));
        uint count = 0;
        Polygonmesh<> m;

        for (uint rid=0; rid<regions.size(); ++rid) {
            std::string base_path = std::string(HOME_PATH) + "../out/" + regions.at(rid) + "/regions_off/";

            // load the region mesh
            std::string region = base_path + "region_" + std::to_string(label) + ".off";
            if (!fs::exists(region)) {
                std::cout << "Warning: region " << region << " does not exist!" << std::endl;
                continue;
            }
            Polygonmesh<> m2(region.c_str());

            // load the region field from file
            std::string cells_data = base_path + "region_" + std::to_string(label) + "_cells_data.csv";
            std::ifstream file(cells_data);
            if (!file.is_open()) {
                std::cerr << "Could not open the file!" << std::endl;
                return 1;
            }
            std::string line;
            std::getline(file, line); // skip the header line
            while (std::getline(file, line)) {
                std::stringstream ss(line);
                std::string cell_id_str, label1_str, label2_str, field_str;
                std::getline(ss, cell_id_str, ',');
                std::getline(ss, label1_str, ',');
                std::getline(ss, label2_str, ',');
                std::getline(ss, field_str, ',');

                uint pid = std::stoi(cell_id_str);
                int l    = std::stod(label1_str);
                double f = std::stod(field_str);
                m2.poly_data(pid).quality = f;
                m2.poly_data(pid).label = l;
            }
            file.close();

            // merge *m2* with *m*, obtaining *res*
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

            std::unordered_map<int, int> tmp_labels_map;
            for (int l=0; l<n_labels; ++l) {
                tmp_labels_map[l] = l;
            }

            std::string output_path = std::string(HOME_PATH) + "out/label_" + std::to_string(label);

            Statistics stats;
            stats.init(output_path, perc_vals, isovals, labels_polys_map, prof);
            stats.compute_stats(res);
            stats.print_stats(res);

            // print verts and cells data (label and field) in a csv file
            print_verts_csv(res, output_path + "/verts_data.csv", prof);
            print_cells_csv(res, tmp_labels_map, output_path + "/cells_data.csv", prof);

            // generate a shapefile containing the iso-regions boundaries
            print_regions_shp(res, stats, labels_polys_map, output_path + "/regions", prof);

            // save a mesh for each iso-region
            print_regions_off(res, tmp_labels_map, n_labels, output_path, prof);

            m = res;
            ++count;
        }
        prof.pop(true, " regions found: " + std::to_string(count));
        std::cout << std::endl;
    }

    return 0;
}
