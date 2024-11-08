#ifndef SEGMENTED_MESH_ATTRIBUTES_H
#define SEGMENTED_MESH_ATTRIBUTES_H

#include <cinolib/meshes/meshes.h>

using namespace cinolib;

struct Vert_std_attributes_double
{
    vec3d          normal  = vec3d(0,0,0);
    Color          color   = Color::WHITE();
    vec3d          uvw     = vec3d(0,0,0);
    int            label   = -1;
    float          quality = 0.0;
    double         fvalue  = 0.0;
    std::bitset<8> flags   = 0x00;
};

struct Poly_std_attributes_double
{
    vec3d          normal  = vec3d(0,0,0);
    Color          color   = Color::WHITE();
    int            label   = -1;
    float          quality = 0.0;
    double         fvalue  = 0.0;
    float          AO      = 1.0;
    std::bitset<8> flags   = 0x00;
};

typedef Mesh_std_attributes M;
typedef Vert_std_attributes_double VD;
typedef Edge_std_attributes E;
typedef Polygon_std_attributes F;
typedef Poly_std_attributes_double PD;

#endif // SEGMENTED_MESH_ATTRIBUTES_H
