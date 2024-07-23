#ifndef CONTOURED_MESH_H
#define CONTOURED_MESH_H

/* class for segmenting a 2D/3D domain. All routines are executed by calling 'run' */
 
#include <cinolib/meshes/meshes.h>
#include <cinolib/color.h>
#include <cinolib/bfs.h>
#include <cinolib/drawable_isocontour.h>
#include <cinolib/drawable_isosurface.h>
#include <cinolib/profiler.h>
#include <cinolib/parallel_for.h>

#include "criteria.h"
#include "statistics.h"
#include "read_parameters.h"

#include <iostream>
#include <fstream>

using namespace cinolib;

class ContouredMesh {

public:
    ContouredMesh(){}
    ~ContouredMesh(){}

    /* complete routine */
    void run(Polygonmesh<> &m, Parameters &Par);
    void run(Polyhedralmesh<> &m, Parameters &Par);

    /* mark edges or faces among cells with different labels */
    template<class M, class V, class E, class P>
    void mark_edges(AbstractMesh<M,E,V,P> &m);
    void mark_faces(Tetmesh<> &m);

    std::vector<double> isovals;

private:
    double field_correction = 0.;
    double field_min, field_max;
    std::map<uint,double> field, field2;
    std::vector<double> field_global;
    std::unordered_map<int, std::vector<uint>> labels_polys_map;
    std::unordered_map<int, int> tmp_labels_map;
    std::unordered_map<int, int> father_son_map;
    Profiler prof;
    bool verbose;

    /* load the scalar field and assign a field value to cells and vertices */
    template<class M, class V, class E, class P> inline
    void load_field(AbstractMesh<M,E,V,P> &m, const std::string field_path);
    void load_global_field(const std::string field_path);
    template<class M, class V, class E, class P> inline
    void load_second_field(AbstractMesh<M,E,V,P> &m, const std::string field_path);

    /* color the mesh elements and vertices wrt the field */
    template<class M, class V, class E, class P> inline
    void color_mesh(AbstractMesh<M,E,V,P> &m);

    /* compute *n_regions+1* isovalues between min_val and max_val */
    void compute_iso(const int n_regions, const std::vector<double> &input_isovals);

    /* insert the isocontours (isosurfaces) cutting the mesh */
    void cut_mesh(Trimesh<> &m);
    void cut_mesh(Tetmesh<> &m);

    /* label the elements inside the isoregions in the mesh */
    template<class M, class V, class E, class P> inline
    void insert_iso1(AbstractMesh<M,E,V,P> &m);
    template<class M, class V, class E, class P> inline
    void insert_iso2(AbstractMesh<M,E,V,P> &m);
    template<class M, class V, class E, class P> inline
    void insert_iso3(AbstractMesh<M,E,V,P> &m);

    /* separate connected components of each label */
    template<class M, class V, class E, class P>
    void separate_ccs(AbstractMesh<M,E,V,P> &m);

    /* removes regions smaller than *FILTER_THRESH* */
    template<class M, class V, class E, class P>
    void remove_small_labels(AbstractMesh<M,E,V,P> &m, const double FILTER_THRESH);

    /* smoothing of the boundaries */
    template<class M, class V, class E, class P>
    uint smooth_boundaries(AbstractMesh<M,E,V,P> &m);

    /* reassign the original *n_regions* labels to the cells */
    template<class M, class V, class E, class P>
    void restore_original_labels(AbstractMesh<M,E,V,P> &m);

    /* create a map <label, all polys with that label> */
    template<class M, class V, class E, class P>
    void update_labels_polys_map(AbstractMesh<M,E,V,P> &m);

    /* save the final mesh and export the cell labels in a csv file */
    template<class M, class V, class E, class P>
    void save_mesh(AbstractMesh<M,E,V,P> &m,
                   const std::string mesh_path, const std::string MESH_FORMAT);
};

#ifndef CINO_STATIC_LIB
#include "contoured_mesh.cpp"
#endif

#endif // CONTOURED_POLYGONMESH_H
