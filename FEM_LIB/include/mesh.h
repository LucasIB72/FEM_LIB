#pragma once
typedef struct Mesh
{
    int dim;

    int n_nodes;
    int n_elements;
    int nodes_per_element;

    int dof_per_node;
    int total_dofs;

    double* coords;        // [n_nodes * dim]

    int* connectivity;     // [n_elements * nodes_per_element]


} Mesh;

void mesh_init(Mesh* mesh);

void mesh_allocate(Mesh* mesh);

void mesh_free(Mesh* mesh);

void mesh_print_info(const Mesh* mesh);