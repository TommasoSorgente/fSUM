#ifndef CONTOURED_MESH_H
#define CONTOURED_MESH_H
 
#include <cinolib/meshes/meshes.h>
#include <cinolib/color.h>
#include <cinolib/bfs.h>
#include <cinolib/drawable_isocontour.h>
#include <cinolib/drawable_isosurface.h>
#include <cinolib/profiler.h>

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

    std::vector<double> isovals;

    template<class M, class E, class V, class P> inline
    void init(AbstractMesh<M,E,V,P> &m, Parameters &Par);

    template<class M, class E, class V, class P> inline
    void segment(AbstractMesh<M,E,V,P> &m, Parameters &Par);

    template<class M, class E, class V, class P> inline
    void filter(AbstractMesh<M,E,V,P> &m, Parameters &Par);

    template<class M, class E, class V, class P> inline
    void smooth_boundaries(AbstractMesh<M,E,V,P> &m, Parameters &Par);

    void output(Polygonmesh<> &m,    Parameters &Par);
    void output(Polyhedralmesh<> &m, Parameters &Par);

private:
    double field_correction = 0.;
    double field_min, field_max;
    std::map<uint,double> field, second_field;
    std::vector<double> global_field;
    std::unordered_map<int, std::vector<uint>> polys_in_region;
    std::unordered_map<int, std::vector<uint>> polys_in_subregion;
    std::unordered_map<int, int> subregion_region_map;
    Profiler prof;
    bool verbose;
    std::string output_path;

    /* load the scalar field and assign a field value to cells and vertices */
    template<class M, class E, class V, class P> inline
    void load_field(AbstractMesh<M,E,V,P> &m, const std::string field_path,
                           std::map<uint,double> &field, bool SMOOTH_FIELD);
    void load_global_field(const std::string field_path);
    template<class M, class E, class V, class P> inline
    void load_second_field(AbstractMesh<M,E,V,P> &m, const std::string field_path);

    /* color the mesh elements and vertices wrt the field */
    template<class M, class E, class V, class P> inline
    void color_mesh(AbstractMesh<M,E,V,P> &m);

    /* compute *n_regions+1* isovalues between min_val and max_val */
    void compute_iso(const int n_regions, const int input_type, const std::vector<double> &input_isovals);

    /* insert the isocontours (isosurfaces) cutting the mesh */
    void cut_mesh(Trimesh<> &m);
    void cut_mesh(Tetmesh<> &m);

    /* label the elements inside the isoregions in the mesh */
    template<class M, class E, class V, class P> inline
    void insert_iso1(AbstractMesh<M,E,V,P> &m);
    template<class M, class E, class V, class P> inline
    void insert_iso2(AbstractMesh<M,E,V,P> &m);
    template<class M, class E, class V, class P> inline
    void insert_iso3(AbstractMesh<M,E,V,P> &m);

    /* separate connected components of each label */
    template<class M, class E, class V, class P> inline
    void separate_ccs(AbstractMesh<M,E,V,P> &m);

    /* removes regions smaller than *FILTER_THRESH* */
    template<class M, class E, class V, class P> inline
    void remove_small_labels(AbstractMesh<M,E,V,P> &m, const double FILTER_THRESH);

    /* smoothing of the boundaries */
    template<class M, class E, class V, class P> inline
    uint smooth_boundaries(AbstractMesh<M,E,V,P> &m);

    /* reassign the original *n_regions* labels to the cells */
    template<class M, class E, class V, class P> inline
    void restore_original_labels(AbstractMesh<M,E,V,P> &m);
};

#ifndef CINO_STATIC_LIB
#include "contoured_mesh.cpp"
#endif

#endif // CONTOURED_POLYGONMESH_H
