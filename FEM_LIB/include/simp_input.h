#ifndef SIMP_INPUT_H
#define SIMP_INPUT_H

#include "mesh.h"

int generate_rect_mesh(
    double L_x, double L_y, double e_s,
    double** coords,
    int** connectivity,
    int* n_nodes,
    int* n_elements,
    int* nodes_per_element);

int generate_rect_mesh_to_mesh(
    double L_x, double L_y, double e_s,
    Mesh* mesh);

#endif
