#include <iostream>
#include <filesystem>

#include <cinolib/meshes/meshes.h>
#include <cinolib/profiler.h>
#include <cinolib/merge_meshes_at_coincident_vertices.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>

#include "../src/segmented_mesh_attributes.h"
#include "../src/statistics.h"
#include "../src/auxiliary.h"
#include "../src/shapefile_utils.h"

using namespace cinolib;
namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
    if (argc < 4) {
        std::cout << "Welcome to Merge Subdomains! Please input:\n"
                     "- the list of subdomains \n"
                     "- the global field \n"
                     "- the isovalues file \n"
                     "- the output directory \n";
        exit(1);
    }

    std::string domains_list   = argv[1];
    std::string field_file     = std::string(HOME_PATH) + argv[2];
    std::string isovalues_file = std::string(HOME_PATH) + argv[3];
    std::string output_dir    = argv[4];

    std::vector<std::string> domains;
    std::stringstream ss(domains_list);
    std::string domain;
    while (ss >> domain) {
        domains.push_back(domain);
    }

    std::vector<double> global_field;
    load_field(field_file, global_field);
    auto minmax = minmax_element(global_field.begin(), global_field.end());
    double field_min = *minmax.first;
    double field_max = *minmax.second;

    std::vector<double> percentiles, isovalues;
    std::vector<int> labels;
    load_isovals(isovalues_file, percentiles, isovalues, labels);

    open_directory(output_dir);

    Profiler prof;
    prof.push("Merge Subdomains");
    GLcanvas *gui = new GLcanvas(1000, 1000);

    for (uint label=0; label<labels.size(); ++label) {
        prof.push("Processed label " + std::to_string(label));
        DrawablePolygonmesh<M,VD,E,PD> *m = new DrawablePolygonmesh<M,VD,E,PD>;
        uint count = 0;
        for (uint dom=0; dom<domains.size(); ++dom) {
            std::string base_path = std::string(HOME_PATH) + "out/segmentation/" + domains.at(dom) + "/regions/";

            // load the domain mesh
            std::string domain = base_path + "region_" + std::to_string(label) + ".obj";
            if (!fs::exists(domain)) {
                std::cout << "Domain " << domain << " does not contain label " << label << std::endl;
                continue;
            }
            DrawablePolygonmesh<M,VD,E,PD> *m_local = new DrawablePolygonmesh<M,VD,E,PD>(domain.c_str());
            std::string cells_data = base_path + "region_" + std::to_string(label) + "_cells_data.csv";
            std::vector<uint> second_label;
            load_cells_data(cells_data, *m_local, second_label); // load field value (fvalue) and label of each cell

            // merge *m* with *m_local*, obtaining *res*
            DrawablePolygonmesh<M,VD,E,PD> *m_merge = new DrawablePolygonmesh<M,VD,E,PD>;
            merge_meshes_at_coincident_vertices(*m, *m_local, *m_merge);
            assert(m_merge->num_polys() == m->num_polys() + m_local->num_polys());
            for (uint pid=0; pid<m->num_polys(); ++pid) {
                m_merge->poly_data(pid) = m->poly_data(pid);
            }
            int poly_offset = m->num_polys();
            for (uint pid=0; pid<m_local->num_polys(); ++pid) {
                m_merge->poly_data(poly_offset + pid) = m_local->poly_data(pid);
            }
            // update m
            m = m_merge;
            ++count;
        }
        prof.pop(true, " regions found: " + std::to_string(count));

        std::unordered_map<int, std::vector<uint>> polys_in_region;
        update_regions_map(*m, polys_in_region);
        std::unordered_map<int, int> subregion_region_map;
        for (int l=0; l<labels.size(); ++l) {
            subregion_region_map[l] = l;
        }

        std::string output_path = output_dir + "/label_" + std::to_string(label);
        open_directory(output_path);
        m->save((output_path + "/mesh.obj").c_str());

        // print verts and cells data (label and field) in a csv file
        print_verts_csv(*m, output_path + "/verts_data.csv", prof);
        print_cells_csv(*m, subregion_region_map, output_path + "/cells_data.csv", prof);

        // print general statistics
        Statistics stats;
        stats.init(percentiles, isovalues, prof);
        stats.compute_stats(*m, polys_in_region);
        stats.print_stats(*m, output_path);

        // generate a shapefile containing the iso-regions boundaries
        print_regions_shp(*m, stats, polys_in_region, output_path + "/shapefile", prof);

        m->mesh_data().filename = "label_" + std::to_string(label);
        for(uint pid=0; pid<m->num_polys(); ++pid) {
            assert(m->poly_data(pid).label == (int)label);
            // float c = (m->poly_data(pid).fvalue - field_min) / (field_max - field_min);
            // float c = (float)m->poly_data(pid).label / (labels.size() - 1);
            // m->poly_data(pid).color = Color::red_white_blue_ramp_01(1. - c);
            m->poly_data(pid).color = Color::parula_ramp(labels.size(), m->poly_data(pid).label);
        }
        // m->poly_color_wrt_label(false);
        m->show_wireframe_transparency(0.1f);
        m->show_marked_edge_color(Color::BLACK());
        m->show_marked_edge_width(2.);
        m->updateGL();
        gui->push(m);
        SurfaceMeshControls<DrawablePolygonmesh<M,VD,E,PD>> *menu = new SurfaceMeshControls<DrawablePolygonmesh<M,VD,E,PD>>(m, gui);
        gui->push(menu);

        std::cout << std::endl;
    }

    prof.pop();
    gui->launch();
    return 0;
}
