#include "mesh.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

void mesh_init(Mesh* mesh)
{
    mesh->dim = 0;

    mesh->n_nodes = 0;
    mesh->n_elements = 0;
    mesh->nodes_per_element = 0;

    mesh->dof_per_node = 0;
    mesh->total_dofs = 0;

    mesh->coords = NULL;
    mesh->connectivity = NULL;
}

void mesh_allocate(Mesh* mesh)
{

    mesh->coords = malloc(
        mesh->n_nodes * mesh->dim * sizeof(double)
    );

    mesh->connectivity = malloc(
        mesh->n_elements *
        mesh->nodes_per_element *
        sizeof(int)
    );

    if (mesh->coords == NULL || mesh->connectivity == NULL)
    {
        printf("Mesh allocation failed\n");
        exit(1);
    }

    mesh->total_dofs = mesh->n_nodes * mesh->dof_per_node;
}

void mesh_free(Mesh* mesh)
{
    if (mesh->coords != NULL)
        free(mesh->coords);

    if (mesh->connectivity != NULL)
        free(mesh->connectivity);

    mesh->coords = NULL;
    mesh->connectivity = NULL;
}

void mesh_print_info(const Mesh* mesh)
{
    printf("Mesh information\n");
    printf("----------------\n");

    printf("Dimension: %d\n", mesh->dim);

    printf("Nodes: %d\n", mesh->n_nodes);

    printf("Elements: %d\n", mesh->n_elements);

    printf("Nodes per element: %d\n", mesh->nodes_per_element);

    printf("DOF per node: %d\n", mesh->dof_per_node);

    printf("Total DOFs: %d\n", mesh->total_dofs);
}