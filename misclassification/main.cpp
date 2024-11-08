#include <iostream>

#include <cinolib/meshes/meshes.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>

#include "../src/segmented_mesh_attributes.h"
#include "../src/auxiliary.h"

using namespace cinolib;

template<class M, class VD, class E, class PD> inline
void load_cells_data(const std::string filepath, Polygonmesh<M,VD,E,PD> &m) {
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
        m.poly_data(pid).fvalue = f;
        m.poly_data(pid).label = l;
    }
    file.close();
}

void load_isovals(const std::string path, std::vector<double> &percentiles, std::vector<double> &isovalues, std::vector<int> &labels) {
    std::ifstream fp(path);
    assert(fp.is_open());
    std::string line;
    std::getline(fp, line);
    std::string p_str, v_str, l_str;
    while (std::getline(fp, line)) {
        std::stringstream ss(line);
        if (std::getline(ss, p_str, ',') && std::getline(ss, v_str, ',')) {
            percentiles.push_back(std::stod(p_str));
            isovalues.push_back(std::stod(v_str));
            if (std::getline(ss, l_str, ',')) {
                labels.push_back(std::stoi(l_str));
            }
        } else {
            std::cerr << "Error reading line: " << line << std::endl;
        }
    }
    fp.close();
}

void load_field(const std::string path, std::vector<double> &field) {
    if (path == "") return;
    std::ifstream fp(path.c_str());
    assert(fp.is_open());
    double f;
    while (fp >> f) {
        field.push_back(f);
    }
    fp.close();
}

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
    if (argc < 3) {
        std::cout << "Welcome to Misclassification Viewer! Please input:\n"
                     "- the mesh \n"
                     "- the mesh cell data \n"
                     "- the isovalues file \n"
                     "and optionally: \n"
                     "- the field_minus_sigma \n"
                     "- the field_plus_sigma \n";
        exit(1);
    }

    std::string mesh_file      = std::string(HOME_PATH) + argv[1];
    std::string cells_file     = std::string(HOME_PATH) + argv[2];
    std::string isovalues_file = std::string(HOME_PATH) + argv[3];
    std::string sigma_m_file = "";
    std::string sigma_p_file = "";
    if (argc == 6) {
        sigma_m_file = std::string(HOME_PATH) + argv[4];
        sigma_p_file = std::string(HOME_PATH) + argv[5];
    }

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

    std::string local_file = std::string(HOME_PATH) + "local_misclassification.txt";
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

    std::string global_file = std::string(HOME_PATH) + "global_misclassification.txt";
    print_global_misclass(m.num_polys(), misclass_vec, n_misclass_vec, global_file);

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
