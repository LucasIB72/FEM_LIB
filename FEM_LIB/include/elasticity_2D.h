#ifndef ELASTICITY_2D_H
#define ELASTICITY_2D_H

#include "element_routine.h"
#include "mesh.h"

typedef struct Material
{
    double E;
    double nu;
} Material;

void elasticity_2d_integrate(double* Ke, double* Re, const double* dN_global, double detJ, double w, const double* material_properties, int n_nodes, int n_dim, int ndof, const double* u_e);

void element_routine_e2d(double* Ke, double* Re, Mesh* mesh, ElementType* etype, Material* mat, int e, const double* u);

#endif
