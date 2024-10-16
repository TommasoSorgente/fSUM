#include <cinolib/gl/glcanvas.h>
#include <cinolib/drawable_isocontour.h>
#include <cinolib/drawable_isosurface.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/gl/volume_mesh_controls.h>
#include <cinolib/meshes/mesh_attributes.h>

#include "src/contoured_mesh.h"
#include "src/read_parameters.h"

using namespace cinolib;

int main(int argc, char *argv[]) {
    string pars_file = (argc == 2) ? argv[1] : string(HOME_PATH) + "parameters.run";
    Read_Parameters Par(pars_file);
    Par.read_file();
    if (argc == 3) { // get mesh and fields from automatic Optimization routine
        Par.set_MESH_PATH(argv[1]);
        Par.set_CLEAN_THRESH(atof(argv[2]));
        Par.set_GUI(false);
    }
    if (argc == 4) { // get mesh and fields from automatic Liguria routine
        Par.set_MESH_PATH(argv[1]);
        Par.set_FIELD_PATH(argv[2]);
        Par.set_FGLOBAL(true);
        Par.set_FIELDG_PATH(argv[3]);
        Par.set_GUI(false);
    }
    if (Par.get_VERBOSE())
        Par.print_pars();

    ContouredMesh FESA;

    switch (Par.get_DIM()) {
    case 2: {
        DrawablePolygonmesh<M,VD,E,PD> m;
        FESA.init(m, Par);
        FESA.segment(m, Par);
        if (Par.get_CLEAN())   FESA.clean(m, Par);
        if (Par.get_SMOOTH())  FESA.smooth(m, Par);
        if (Par.get_ANALYZE()) FESA.output(m, Par);
        else                   FESA.restore_original_labels(m);
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
            std::vector<double> isovalues = FESA.isovals;
            for (uint i=0; i<isovalues.size(); ++i) {
                DrawableIsocontour<M,VD,E,PD> *contour = new DrawableIsocontour<M,VD,E,PD>(m, isovalues.at(i));
                contour->thickness = 30000; //700;
                contour->color = Color::BLACK(); //parula_ramp(isovalues.size(), i);
                gui.push(contour);
            }
        }

        bool MISCLASS = false;
        if (MISCLASS) {
            std::string misclass_path = FESA.output_path + "/misclassifications.csv";
            std::ifstream fp(misclass_path.c_str());
            assert(fp.is_open());

            std::string line;
            std::getline(fp, line);

            uint pid;
            vec3d c;
            double f, l_min, l_max, s_min, s_max, misclass;
            char d;
            while (std::getline(fp, line)) {
                std::stringstream ss(line);
                ss >> pid >> d >> c.x() >> d >> c.y() >> d >> c.z() >> d >> f >>
                      d >> l_min >> d >> l_max >> d >> s_min >> d >> s_max >> d >> misclass;
                gui.push_marker(c, "", Color::BLACK(), 4.);
            }
            fp.close();
        }

        SurfaceMeshControls<DrawablePolygonmesh<M,VD,E,PD>> menu(&m, &gui);
        gui.push(&menu);
        gui.camera.rotate_x(-90);
        gui.update_GL_matrices();
        gui.launch();
        break;
    }
    case 3: {
        DrawablePolyhedralmesh<M,VD,E,F,PD> m;
        FESA.init(m, Par);
        FESA.segment(m, Par);
        if (Par.get_CLEAN())   FESA.clean(m, Par);
        if (Par.get_SMOOTH())  FESA.smooth(m, Par);
        if (Par.get_ANALYZE()) FESA.output(m, Par);
        else                   FESA.restore_original_labels(m);
        if (!Par.get_GUI()) break;

        m.show_out_wireframe_transparency(0.2f);
        int n_labels = Par.get_N_REGIONS() - 1;
        for(uint pid=0; pid<m.num_polys(); ++pid) {
            float c = (float)m.poly_data(pid).label / n_labels;
            m.poly_data(pid).color = Color::red_white_blue_ramp_01(1. - c);
        }
        // m.poly_color_wrt_label(true);
        m.show_out_poly_color();
        mark_faces(m);
        m.updateGL();

        GLcanvas gui(1000, 1000);
        gui.push(&m);

        if (Par.get_ISOCONTOURS()) assert(false && "show_iso not working for 3D meshes!");

        VolumeMeshControls<DrawablePolyhedralmesh<M,VD,E,F,PD>> menu(&m, &gui);
        gui.push(&menu);
        gui.launch();
        break;
    }
    default:
        cout << "ERROR: wrong dimension!\n";
    }
}
