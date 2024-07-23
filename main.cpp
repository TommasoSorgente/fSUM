#include <cinolib/gl/glcanvas.h>
#include <cinolib/drawable_isocontour.h>
#include <cinolib/drawable_isosurface.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/gl/volume_mesh_controls.h>

#include "src/contoured_mesh.h"
#include "src/read_parameters.h"

using namespace cinolib;

int main(int argc, char *argv[]) {
    string pars_file = (argc == 2) ? argv[1] : string(HOME_PATH) + "parameters.run";
    Read_Parameters Par(pars_file);
    Par.read_file();
    if (argc == 4) { // get mesh and fields from automatic routine
        Par.set_MESH_PATH(argv[1]);
        Par.set_FIELD1_PATH(argv[2]);
        Par.set_FIELD2_PATH(argv[3]);
        Par.set_GUI(false);
    }
    if (Par.get_VERBOSE())
        Par.print_pars();

    ContouredMesh Cont_mesh;

    switch (Par.get_DIM()) {
    case 2: {
        DrawablePolygonmesh<> m;
        Cont_mesh.run(m, Par);
        if (!Par.get_GUI()) break;

        m.show_wireframe_transparency(0.2f);
        int n_labels = Par.get_N_REGIONS() - 1;
        for(uint pid=0; pid<m.num_polys(); ++pid) {
            float c = (float)m.poly_data(pid).label / n_labels;
            m.poly_data(pid).color = Color::red_white_blue_ramp_01(1. - c);
        }
        // m.poly_color_wrt_label(true);
        m.show_poly_color();
        // m.show_vert_color();
        m.updateGL();

        GLcanvas gui(1000, 1000);
        gui.push(&m);

        if (Par.get_SHOW_ISO()) {
            std::vector<double> isovalues = Cont_mesh.isovals;
            for (uint i=0; i<isovalues.size(); ++i) {
                DrawableIsocontour<> *contour = new DrawableIsocontour<>(m, isovalues.at(i));
                contour->thickness = 700;
                contour->color = Color::BLACK(); //parula_ramp(isovalues.size(), i);
                gui.push(contour);
            }
        }

        SurfaceMeshControls<DrawablePolygonmesh<>> menu(&m, &gui);
        gui.push(&menu);
        // gui.camera.rotate_x(-90);
        // gui.update_GL_matrices();
        gui.launch();
        break;
    }
    case 3: {
        DrawablePolyhedralmesh<> m;
        Cont_mesh.run(m, Par);
        if (!Par.get_GUI()) break;

        m.show_out_wireframe_transparency(0.2f);
        m.poly_color_wrt_label(true);
        m.show_out_poly_color();
        m.updateGL();

        GLcanvas gui(1000, 1000);
        gui.push(&m);

        if (Par.get_SHOW_ISO()) {
            assert(false && "show_iso not working for 3D meshes!");
            std::vector<double> isovalues = Cont_mesh.isovals;
            Tetmesh<> m_tet(m.vector_verts(), m.vector_polys());
            for (uint i=0; i<isovalues.size(); ++i) {
                DrawableIsosurface<> *surface = new DrawableIsosurface<>(m_tet, isovalues.at(i));
                surface->color = Color::BLACK(); //parula_ramp(isovalues.size(), i);
                gui.push(surface);
            }
        }

        VolumeMeshControls<DrawablePolyhedralmesh<>> menu(&m, &gui);
        gui.push(&menu);
        gui.launch();
        break;
    }
    default:
        cout << "ERROR: wrong dimension!\n";
    }
}
