#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "mesh.h"
#include "elements.h"

void assemble_global_stiffness(
    double* K_global,
    double* R_global,
    Mesh* mesh,
    ElementType* etype,
    PhysicsModel* physics,
    const double* material_properties,
    const double* u);

void assemble_global_stiffness_sparse(
    int** rowIndex,
    int** columns,
    double** values,
    int* nnz,
    Mesh* mesh,
    ElementType* etype,
    PhysicsModel* physics,
    const double* material_properties);

#endif
