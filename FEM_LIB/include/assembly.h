#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "mesh.h"
#include "elements.h"
#include "elasticity_2D.h"

void assemble_global_stiffness(
    double* K_global,
    double* R_global,
    Mesh* mesh,
    ElementType* etype,
    Material* mat,
    const double* u);

void assemble_global_stiffness_sparse(
    int** rowIndex,
    int** columns,
    double** values,
    int* nnz,
    Mesh* mesh,
    ElementType* etype,
    Material* mat);

#endif
