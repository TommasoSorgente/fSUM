#include <cinolib/gl/glcanvas.h>
#include <cinolib/drawable_isocontour.h>
#include <cinolib/drawable_isosurface.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/gl/volume_mesh_controls.h>
#include <cinolib/meshes/mesh_attributes.h>
#include <tclap/CmdLine.h>

#include "src/segmented_mesh.h"
#include "src/parameters.h"

using namespace cinolib;

int main(int argc, char *argv[]) {
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
        else                   S.restore_original_labels(m);
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
