#include <cinolib/gl/glcanvas.h>
#include <cinolib/drawable_isocontour.h>
#include <cinolib/drawable_isosurface.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/gl/volume_mesh_controls.h>
#include <cinolib/meshes/mesh_attributes.h>
#include <tclap/CmdLine.h>

#include "../src/segmented_mesh.h"
#include "../src/parameters.h"

using namespace cinolib;

int main(int argc, char *argv[]) {
    std::string data_folder = "../../data/Liguria_grid/";
    std::vector<std::string> domains = {"1_5TERRE", "2_MAGRA", "3_PETRONIO", "4_PORTOFINO", "5_ENTELLA", "6_AVETO", "7_PADANO", "8_BISAGNO", "9_POLCEVERA", "10_ARENZANO", "11_SASSELLO", "12_BORMIDE", "13_SAVONESE", "14_IMPERIESE"};

    GLcanvas gui(1000, 1000);
    std::vector<double> global_field;
    double field_min, field_max;

    for (std::string s : domains) {
        std::string mesh_path = data_folder + s + "/mesh.obj";
        DrawablePolygonmesh<M,VD,E,PD> *m = new DrawablePolygonmesh<M,VD,E,PD>(mesh_path.c_str());

        std::string field_path = data_folder + s + "/field.csv";
        std::map<uint,double> field;
        std::ifstream fp(field_path.c_str());
        assert(fp.is_open());
        double f;
        int i = 0;
        while (fp >> f) {
            field[i] = f;
            ++i;
        }
        fp.close();
        assert(field.size() == m->num_polys());
        for (uint pid=0; pid<m->num_polys(); ++pid) {
            m->poly_data(pid).fvalue = field.at(pid);
        }

        for (uint i=0; i<field.size(); ++i) {
            global_field.push_back(field[i]);
        }
        auto minmax = minmax_element(global_field.begin(), global_field.end());
        field_min = *minmax.first;
        field_max = *minmax.second;

        // set vert_data.uvw.u from the weighted average of poly.data.fvalue
        for (uint vid=0; vid<m->num_verts(); ++vid) {
            double value = 0., mass  = 0.;
            for (uint pid : m->adj_v2p(vid)) {
                double coeff = m->poly_mass(pid) / m->verts_per_poly(pid);
                value += field.at(pid) * coeff;
                mass  += coeff;
            }
            if (mass < 1e-4)
                m->vert_data(vid).fvalue = field.at(m->adj_v2p(vid).front());
            // assert(mass > 0);
            else
                m->vert_data(vid).fvalue = value / mass;
        }

        for (uint pid=0; pid<m->num_polys(); ++pid) {
            // double c = (field.at(pid) - field_min) / (field_max - field_min); // input field
            double c = (m->poly_data(pid).fvalue - field_min) / (field_max - field_min); // mean field
            m->poly_data(pid).color = Color::red_white_blue_ramp_01(1. - c);
        }

        for (uint vid=0; vid<m->num_verts(); ++vid) {
            double c = (m->vert_data(vid).fvalue - field_min) / (field_max - field_min);
            m->vert_data(vid).color = Color::red_white_blue_ramp_01(1. - c);
            m->vert_data(vid).uvw[0] = m->vert_data(vid).fvalue; // for drawing isocontours
        }

        m->show_wireframe_transparency(0.1f);
        m->show_poly_color();
        m->show_marked_edge_width(2.);
        m->show_marked_edge_color(Color::GREEN());
        m->updateGL();

        SurfaceMeshControls<DrawablePolygonmesh<M,VD,E,PD>> *menu = new SurfaceMeshControls<DrawablePolygonmesh<M,VD,E,PD>>(m, &gui);
        gui.push(menu);
        gui.push(m);
    }
    gui.launch();
}

int main2(int argc, char *argv[]) {
    Parameters Par;
    try {
        TCLAP::CmdLine cmd("S", ' ', "0.9");
        Par.read_cmdline(cmd, argc, argv);
    }
    catch (exception &e) {
        std::cerr << e.what() << std::endl;
        exit(1);
    }
    if (Par.get_VERBOSE()) Par.print_pars();

    SegmentedMesh S;

    switch (Par.get_DIM()) {
    case 2: {
        DrawablePolygonmesh<M,VD,E,PD> m;
        S.init(m, Par);
        S.segment(m, Par);
        if (Par.get_CLEAN())   S.clean(m, Par);
        if (Par.get_SMOOTH())  S.smooth(m, Par);
        if (Par.get_ANALYZE()) S.output(m, Par);
        S.restore_original_labels(m);
        if (!Par.get_GUI()) break;

        m.show_wireframe_transparency(0.1f);
        int n_labels = Par.get_N_REGIONS() - 1;
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

        GLcanvas gui(1000, 1000);
        gui.push(&m);

        if (Par.get_ISOCONTOURS()) {
            std::vector<double> isovalues = S.isovals;
            for (uint i=0; i<isovalues.size(); ++i) {
                DrawableIsocontour<M,VD,E,PD> *contour = new DrawableIsocontour<M,VD,E,PD>(m, isovalues.at(i));
                contour->thickness = 30000; //700;
                contour->color = Color::BLACK(); //parula_ramp(isovalues.size(), i);
                gui.push(contour);
            }
        }

        SurfaceMeshControls<DrawablePolygonmesh<M,VD,E,PD>> menu(&m, &gui);
        gui.push(&menu);
        // gui.camera.rotate_x(-90);
        // gui.update_GL_matrices();
        gui.launch();
        break;
    }
    case 3: {
        DrawablePolyhedralmesh<M,VD,E,F,PD> m;
        S.init(m, Par);
        S.segment(m, Par);
        if (Par.get_CLEAN())   S.clean(m, Par);
        if (Par.get_SMOOTH())  S.smooth(m, Par);
        if (Par.get_ANALYZE()) S.output(m, Par);
        S.restore_original_labels(m);
        if (!Par.get_GUI()) break;

        m.show_out_wireframe_transparency(0.1f);
        int n_labels = Par.get_N_REGIONS() - 1;
        for(uint pid=0; pid<m.num_polys(); ++pid) {
            float c = (float)m.poly_data(pid).label / n_labels;
            // m.poly_data(pid).color = Color::red_white_blue_ramp_01(1. - c);
            m.poly_data(pid).color = Color::parula_ramp(Par.get_N_REGIONS(), m.poly_data(pid).label);
        }
        // m.poly_color_wrt_label(true);
        m.show_out_poly_color();
        mark_faces(m);
        m.updateGL();

        GLcanvas gui(1000, 1000);
        gui.push(&m);

        if (Par.get_ISOCONTOURS()) assert(false && "ISOCONTOURS not available for 3D meshes!");

        VolumeMeshControls<DrawablePolyhedralmesh<M,VD,E,F,PD>> menu(&m, &gui);
        gui.push(&menu);
        gui.launch();
        break;
    }
    default:
        cout << "ERROR: wrong dimension!\n";
    }
}
