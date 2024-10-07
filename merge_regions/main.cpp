#include <iostream>
#include <filesystem>

#include <cinolib/meshes/meshes.h>
#include <cinolib/profiler.h>
#include <cinolib/merge_meshes_at_coincident_vertices.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>

#include "../src/read_parameters.h"
#include "../src/statistics.h"
#include "../src/auxiliary.h"
#include "../src/shapefile_utils.h"

using namespace cinolib;
namespace fs = std::filesystem;

void compute_isovals(const std::vector<double> &global_field, const uint n_labels, const std::vector<double> &input_vals, const uint type, std::vector<double> &isovals, std::vector<double> &percentiles)
{
    auto minmax = minmax_element(global_field.begin(), global_field.end());
    double field_min = *minmax.first;
    double field_max = *minmax.second;

    switch (type) {
    case 1: { // equi-spaced isovals
        double step = abs(field_max - field_min) / n_labels;
        for (uint i = 0; i < n_labels; ++i) {
            double isoval = field_min + step * i;
            isovals.push_back(isoval);
            double perc = inv_percentile(global_field, isoval);
            percentiles.push_back(perc);
        }
        isovals.push_back(field_max);
        break;
    }
    case 2: { // percentile isovals
        assert(input_vals.front()==0 && input_vals.back()==100);
        percentiles = input_vals;
        isovals.push_back(field_min);
        for (int i=1; i<input_vals.size()-1; ++i) {
            int perc = input_vals.at(i);
            double isoval = percentile(global_field, perc);
            isovals.push_back(isoval);
        }
        isovals.push_back(field_max);
        break;
    }
    case 3: { // manually set isovals
        for (uint i = 0; i < input_vals.size()-1; ++i) {
            assert(input_vals.at(i) < input_vals.at(i+1) && "isovalues are not ordered");
        }
        for (double isoval : input_vals) {
            double perc = inv_percentile(global_field, isoval);
            percentiles.push_back(perc);
        }
        isovals = input_vals;
        break;
    }
    default:
        break;
    }
    assert((int)isovals.size() == n_labels+1 && "#isovalues != #regions+1");
}

void load_cells_data(const string filepath, Polygonmesh<> &m) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Could not open the file " << filepath << std::endl;
    }
    std::string line;
    std::getline(file, line); // skip the header line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell_id_str, x, y, z, label1_str, label2_str, field_str;
        std::getline(ss, cell_id_str, ',');
        std::getline(ss, x, ',');
        std::getline(ss, y, ',');
        std::getline(ss, z, ',');
        std::getline(ss, label1_str, ',');
        std::getline(ss, label2_str, ',');
        std::getline(ss, field_str, ',');

        uint pid = std::stoi(cell_id_str);
        int l    = std::stod(label1_str);
        double f = std::stod(field_str);
        m.poly_data(pid).quality = f;
        m.poly_data(pid).label = l;
    }
    file.close();
}

int main(int argc, char *argv[])
{
    Profiler prof;
    prof.push("Merge Regions");

    Read_Parameters Par(string(HOME_PATH) + "../parameters.run");
    Par.read_file();

    std::vector<std::string> domains{
                                     "1_F5TERRE",
                                     "2_FMAGRA",
                                     "3_FPETRONIO",
                                     "4_FPORTOFINO",
                                     "5_FENTELLA",
                                     "6_FAVETO",
                                     "7_FPADANO",
                                     "8_FBISAGNO",
                                     "9_FPOLCEVERA",
                                     "10_FARENZANO",
                                     "11_FSASSELLO",
                                     "12_FBORMIDE",
                                     "13_FSAVONESE",
                                     "14_FIMPERIESE"
                                    };
    std::string field_path = std::string(HOME_PATH) + "../data/Liguria_grid/cr_mean_global.csv";

    std::vector<double> global_field;
    std::ifstream fp(field_path);
    assert(fp.is_open());
    double f;
    while (fp >> f) {
        global_field.push_back(f);
    }
    fp.close();

    uint n_labels = Par.get_N_REGIONS();
    std::vector<double> isovals;
    std::vector<double> percentiles;
    compute_isovals(global_field, n_labels, Par.get_ISOVAL_VALS(), Par.get_ISOVAL_TYPE(), isovals, percentiles);

    GLcanvas *gui = new GLcanvas(1000, 1000);

    for (uint label=0; label<n_labels; ++label) {
        prof.push("Processed label " + std::to_string(label));
        DrawablePolygonmesh<> *m = new DrawablePolygonmesh<>;
        uint count = 0;
        for (uint dom=0; dom<domains.size(); ++dom) {
            std::string base_path = std::string(HOME_PATH) + "../out/" + domains.at(dom) + "/regions/";

            // load the domain mesh
            std::string domain = base_path + "region_" + std::to_string(label) + ".obj";
            if (!fs::exists(domain)) {
                std::cout << "Domain " << domain << " does not contain label " << label << std::endl;
                continue;
            }
            DrawablePolygonmesh<> *m_local = new DrawablePolygonmesh<>(domain.c_str());
            std::string cells_data = base_path + "region_" + std::to_string(label) + "_cells_data.csv";
            load_cells_data(cells_data, *m_local); // load field value (quality) and label of each cell

            // merge *m* with *m_local*, obtaining *res*
            DrawablePolygonmesh<> *m_merge = new DrawablePolygonmesh<>;
            merge_meshes_at_coincident_vertices(*m, *m_local, *m_merge);
            for (uint pid=0; pid<m->num_polys(); ++pid) {
                m_merge->poly_data(pid) = m->poly_data(pid);
            }
            int offset = m->num_polys();
            for (uint pid=0; pid<m_local->num_polys(); ++pid) {
                m_merge->poly_data(offset + pid) = m_local->poly_data(pid);
            }

            // update m
            m = m_merge;
            ++count;
        }
        prof.pop(true, " regions found: " + std::to_string(count));

        std::unordered_map<int, std::vector<uint>> polys_in_region;
        update_regions_map(*m, polys_in_region);
        std::unordered_map<int, int> subregion_region_map;
        for (int l=0; l<n_labels; ++l) {
            subregion_region_map[l] = l;
        }

        std::string output_path = std::string(HOME_PATH) + "out/label_" + std::to_string(label);
        open_directory(output_path);
        m->save((output_path + "/mesh.obj").c_str());

        // print verts and cells data (label and field) in a csv file
        print_verts_csv(*m, output_path + "/verts_data.csv", prof);
        print_cells_csv(*m, subregion_region_map, output_path + "/cells_data.csv", prof);

        // print general statistics
        Statistics stats;
        stats.init(percentiles, isovals, prof);
        stats.compute_stats(*m, polys_in_region);
        stats.print_stats(*m, output_path);

        // generate a shapefile containing the iso-regions boundaries
        print_regions_shp(*m, stats, polys_in_region, output_path + "/shapefile", prof);

        if (Par.get_GUI()) {
            m->mesh_data().filename = "label_" + std::to_string(label);
            int n_labels = Par.get_N_REGIONS() - 1;
            for(uint pid=0; pid<m->num_polys(); ++pid) {
                float c = (float)m->poly_data(pid).label / n_labels;
                m->poly_data(pid).color = Color::red_white_blue_ramp_01(1. - c);
            }
            // m->poly_color_wrt_label(false);
            m->updateGL();
            gui->push(m);
            SurfaceMeshControls<DrawablePolygonmesh<>> *menu = new SurfaceMeshControls<DrawablePolygonmesh<>>(m, gui);
            gui->push(menu);
        }
        std::cout << std::endl;
    }

    prof.pop();
    if (Par.get_GUI())
        gui->launch();
    return 0;
}
