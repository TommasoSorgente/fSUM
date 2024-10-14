#ifndef CRITERIA_H
#define CRITERIA_H

/* class containing the criteria used by 'contoured_mesh' for classifying regions */
 
#include <cinolib/meshes/meshes.h>

using namespace cinolib;

/* is the variance of the field among the cells in *polys* smaller than threshold? */
template<class M, class V, class E, class P> inline
double CRIT_variance(const AbstractMesh<M,V,E,P> &m,
                     const std::vector<uint> &polys, const double THRESH)
{
    double sum = 0.;
    std::unordered_set<double> field;
    for (uint pid : polys) {
        double f = m.poly_data(pid).fvalue;
        field.insert(f);
        sum += f;
    }
    double mean = sum / polys.size();
    double var = 0.;
    for(auto f : field) {
      var += (f - mean) * (f - mean);
    }
    var = sqrt(var / polys.size());
    assert(0. <= var && var <= 1.);
    return (var < THRESH);
}

/**********************************************************************/

/* is the mass (area or volume) covered by the cells in *polys* smaller than threshold? */
template<class M, class V, class E, class P> inline
bool CRIT_mass(const AbstractMesh<M,V,E,P> &m, const std::vector<uint> &polys, const double THRESH)
{
    double mass = 0.;
    for (uint pid : polys) {
        mass += m.poly_mass(pid);
        if (mass > THRESH)
            return false;
    }
    return (mass < THRESH);
}

/**********************************************************************/

/* is the shape ratio of the mesh smaller than threshold?
 * (shape ratio is defined to reach 1 on the circle and the sphere) */
template<class M, class V, class E, class P> inline
bool CRIT_shape(const Polygonmesh<M,V,E,P> &mesh, const double THRESH)
{
    double area = mesh.mesh_area();
    double perimeter = 0.;
    for (auto e : mesh.get_boundary_edges()) {
        vec3d v0 = mesh.vert(e.first);
        vec3d v1 = mesh.vert(e.second);
        perimeter += v0.dist(v1);
    }
    double shape_ratio = (4. * M_PI * area) / (perimeter * perimeter);
    assert(0. <= shape_ratio && shape_ratio <= 1.);
    return (shape_ratio < THRESH);
}

template<class M, class V, class E, class P> inline
bool CRIT_shape(const Trimesh<M,V,E,P> &_mesh, const double THRESH)
{
    Polygonmesh<M,V,E,P> mesh(_mesh.vector_verts(), _mesh.vector_polys());
    return CRIT_shape(mesh, THRESH);
}

template<class M, class V, class E, class P> inline
bool CRIT_shape(const Polyhedralmesh<M,V,E,P> &mesh, const double THRESH)
{
    double volume = mesh.mesh_volume();
    double srf_area = mesh.mesh_srf_area();
    double shape_ratio = (36. * M_PI * volume * volume) / (srf_area * srf_area * srf_area);
    assert(0. <= shape_ratio && shape_ratio <= 1.);
    return (shape_ratio < THRESH);
}

template<class M, class V, class E, class P> inline
bool CRIT_shape(const Tetmesh<M,V,E,P> &_mesh, const double THRESH)
{
    std::vector<std::vector<bool>> windings(_mesh.num_polys());
    for (uint pid = 0; pid < _mesh.num_polys(); ++pid) {
        std::vector<bool> w(_mesh.faces_per_poly(pid), true);
        windings.at(pid) = w;
    }
    Polyhedralmesh<M,V,E,P> mesh(_mesh.vector_verts(), _mesh.vector_faces(), _mesh.vector_polys(), windings);
    return CRIT_shape(mesh, THRESH);
}

#endif // CRITERIA_H
