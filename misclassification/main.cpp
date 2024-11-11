#include <iostream>

#include <cinolib/meshes/meshes.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>

#include "../src/segmented_mesh_attributes.h"
#include "../src/auxiliary.h"

using namespace cinolib;

void print_global_misclass(const uint n_polys, const std::vector<double> misclass_vec, const std::vector<uint> n_misclass_vec, const std::string path) {
    std::ofstream fp(path.c_str());
    assert(fp.is_open());
    fp << "# n regions, misclassification, % misclassified, avg misclassification" << std::endl;

    double regions = misclass_vec.size();
    double misclass = std::accumulate(misclass_vec.begin(), misclass_vec.end(), 0.);
    double n = std::accumulate(n_misclass_vec.begin(), n_misclass_vec.end(), 0.);
    double avg_misclass = misclass / n;
    n /= n_polys;

    fp << regions << ", " << misclass << ", " << n << ", " << avg_misclass << std::endl;
    fp.close();
}

int main(int argc, char *argv[])
{
    if (argc < 5) {
        std::cout << "Welcome to Misclassification Viewer! Please input:\n"
                     "- the mesh \n"
                     "- the mesh cell data \n"
                     "- the isovalues file \n"
                     "- the output path \n"
                     "- the gui (1: yes, 0: no) \n"
                     "and optionally: \n"
                     "- the field_minus_sigma \n"
                     "- the field_plus_sigma \n";
        exit(1);
    }

    std::string mesh_file      = std::string(HOME_PATH) + argv[1];
    std::string cells_file     = std::string(HOME_PATH) + argv[2];
    std::string isovalues_file = std::string(HOME_PATH) + argv[3];
    std::string output_path    = argv[4];
    bool use_gui               = atoi(argv[5]);
    std::string sigma_m_file = "";
    std::string sigma_p_file = "";
    if (argc == 8) {
        sigma_m_file = std::string(HOME_PATH) + argv[6];
        sigma_p_file = std::string(HOME_PATH) + argv[7];
    }
    open_directory(output_path);

    DrawablePolygonmesh<M,VD,E,PD> m(mesh_file.c_str());
    load_cells_data(cells_file, m);

    std::vector<double> percentiles, isovalues;
    std::vector<int> labels;
    load_isovals(isovalues_file, percentiles, isovalues, labels);

    std::vector<double> field_m_sigma, field_p_sigma;
    load_field(sigma_m_file, field_m_sigma);
    load_field(sigma_p_file, field_p_sigma);

    std::unordered_map<int, std::vector<uint>> polys_in_region;
    update_regions_map(m, polys_in_region);

    std::vector<double> misclass_vec;
    std::vector<uint> n_misclass_vec(1,0);
    std::vector<uint> misclass_pids;

    std::string local_file = output_path + "/local_misclassification.txt";
    std::ofstream fp(local_file.c_str());
    assert(fp.is_open());
    fp << "# pid, centroid x, centroid y, centroid z, field, isovalue_min, isovalue_max, field_m_sigma, field_p_sigma, misclassification" << std::endl;
    fp << std::fixed;

    for (auto &l : polys_in_region) {
        double misclass = 0.;

        uint lid = *std::find(labels.begin(), labels.end(), l.first);
        std::pair<double, double> lambda(isovalues.at(lid), isovalues.at(lid+1));

        for (uint pid : l.second) {
            std::pair<double,double> sigma = {m.poly_data(pid).fvalue, m.poly_data(pid).fvalue};
            if (field_m_sigma.size() > 0 && field_p_sigma.size() > 0) {
                sigma = {field_m_sigma[pid], field_p_sigma[pid]};
            }
            double d = poly_misclassification(m, pid, lambda, sigma);
            vec3d c = m.poly_centroid(pid);
            if (d > 1e-4) {
                fp.precision(2);
                fp << pid << ", " << c.x() << ", " << c.y() << ", " << c.z() << ", ";
                fp.precision(6);
                fp << m.poly_data(pid).fvalue << ", "
                   << lambda.first << ", " << lambda.second << ", "
                   << sigma.first << ", " << sigma.second << ", "
                   << d << std::endl;
                misclass += d;
                ++n_misclass_vec.back();
                misclass_pids.push_back(pid);
            }
        }
        misclass_vec.push_back(misclass);
    }
    fp.close();

    std::string global_file = output_path + "/global_misclassification.txt";
    print_global_misclass(m.num_polys(), misclass_vec, n_misclass_vec, global_file);

    if (!use_gui)
        exit(0);

    GLcanvas gui(1000, 1000);
    gui.push(&m);
    m.show_wireframe_transparency(0.1f);
    int n_labels = labels.size() - 1;
    for(uint pid=0; pid<m.num_polys(); ++pid) {
        float c = (float)m.poly_data(pid).label / n_labels;
        m.poly_data(pid).color = Color::red_white_blue_ramp_01(1. - c);
    }
    // m.poly_color_wrt_label(false);
    m.show_poly_color();
    mark_edges(m);
    m.show_marked_edge_width(1.);
    m.show_marked_edge_color(Color::BLACK());
    m.updateGL();

    for (uint pid : misclass_pids) {
        gui.push_marker(m.poly_centroid(pid), "", Color::BLACK(), 3.);
    }
    fp.close();

    SurfaceMeshControls<DrawablePolygonmesh<M,VD,E,PD>> menu(&m, &gui);
    gui.push(&menu);
    gui.launch();
}
