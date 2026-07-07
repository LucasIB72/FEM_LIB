#ifndef ELASTICITY_2D_PLANE_STRESS_H
#define ELASTICITY_2D_PLANE_STRESS_H

#include "element_routine.h"
#include "mesh.h"

typedef struct Material
{
    double E;
    double nu;
} Material;

void elasticity_2d_plane_stress_integrate(double* Ke, double* Re, const double* dN_global, double detJ, double w, const double* material_properties, int n_nodes, int n_dim, int ndof, const double* u_e);

#endif
