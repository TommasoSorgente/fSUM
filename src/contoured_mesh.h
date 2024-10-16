#ifndef CONTOURED_MESH_H
#define CONTOURED_MESH_H
 
#include <cinolib/meshes/meshes.h>
#include <cinolib/color.h>
#include <cinolib/bfs.h>
#include <cinolib/drawable_isocontour.h>
#include <cinolib/drawable_isosurface.h>
#include <cinolib/profiler.h>

#include "contoured_mesh_attributes.h"
#include "criteria.h"
#include "statistics.h"
#include "auxiliary.h"
#include "shapefile_utils.h"
#include "read_parameters.h"

#include <iostream>
#include <fstream>

using namespace cinolib;

/* class for segmenting a 2D/3D domain */
class ContouredMesh {

public:
    ContouredMesh(){}
    ~ContouredMesh(){}

    std::vector<double> percentiles, isovals;

    template<class M, class V,  class E, class P> inline
    void init(AbstractMesh<M,V,E,P> &m, Parameters &Par);

    template<class M, class V,  class E, class P> inline
    void segment(AbstractMesh<M,V,E,P> &m, Parameters &Par);

    template<class M, class V,  class E, class P> inline
    void clean(AbstractMesh<M,V,E,P> &m, Parameters &Par);

    template<class M, class V,  class E, class P> inline
    void smooth(AbstractMesh<M,V,E,P> &m, Parameters &Par);

    template<class M, class V,  class E, class P> inline
    void output(Polygonmesh<M,V,E,P> &m,    Parameters &Par);
    template<class M, class E, class V, class F, class P> inline
    void output(Polyhedralmesh<M,E,V,F,P> &m, Parameters &Par);
    std::string output_path;

    /* reassign the original *n_regions* labels to the cells */
    template<class M, class V,  class E, class P> inline
    void restore_original_labels(AbstractMesh<M,V,E,P> &m);

private:
    double field_correction = 0.;
    double field_min, field_max;
    std::map<uint,double> field;
    std::vector<double> global_field;
    std::unordered_map<int, std::vector<uint>> polys_in_region;
    std::unordered_map<int, std::vector<uint>> polys_in_subregion;
    std::unordered_map<int, int> subregion_region_map;
    Profiler prof;
    bool verbose;

    /* load the scalar field and assign a field value to cells and vertices */
    template<class M, class V,  class E, class P> inline
    void load_field(AbstractMesh<M,V,E,P> &m, const std::string field_path, const int field_type);
    template<class M, class V,  class E, class P> inline
    void extend_field_to_verts(AbstractMesh<M,V,E,P> &m);
    template<class M, class V,  class E, class P> inline
    void extend_field_to_cells(AbstractMesh<M,V,E,P> &m);

    void load_global_field(const std::string field_path);

    /* color the mesh elements and vertices wrt the field */
    template<class M, class V,  class E, class P> inline
    void color_mesh(AbstractMesh<M,V,E,P> &m);

    /* compute *n_regions+1* isovalues between min_val and max_val */
    void compute_isovalues(const int n_regions, const int input_type, const std::vector<double> &input_isovals);

    /* insert the isocontours (isosurfaces) cutting the mesh */
    void cut_mesh(Trimesh<> &m);
    void cut_mesh(Tetmesh<> &m);

    /* label the elements inside the isoregions in the mesh */
    template<class M, class V,  class E, class P> inline
    void compute_isoregions(AbstractMesh<M,V,E,P> &m, const bool DENOISE);
    template<class M, class V,  class E, class P> inline
    void compute_isoregions2(AbstractMesh<M,V,E,P> &m, const bool DENOISE);

    /* separate connected components of each label */
    template<class M, class V,  class E, class P> inline
    void separate_ccs(AbstractMesh<M,V,E,P> &m);

    /* removes regions smaller than *FILTER_THRESH* */
    template<class M, class V,  class E, class P> inline
    void remove_small_labels(AbstractMesh<M,V,E,P> &m, const double FILTER_THRESH);

    /* smoothing of the boundaries */
    template<class M, class V,  class E, class P> inline
    uint smooth_boundaries(AbstractMesh<M,V,E,P> &m);
};

#ifndef CINO_STATIC_LIB
#include "contoured_mesh.cpp"
#endif

#endif // CONTOURED_POLYGONMESH_H
