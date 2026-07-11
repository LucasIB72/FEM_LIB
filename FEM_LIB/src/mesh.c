#include "mesh.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

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

    mesh->n_groups = 0;
    mesh->groups = NULL;
}

void mesh_allocate(Mesh* mesh)
{
    //Aloca a memoria do vetor de coordenadas (stride 3: x, y, z mesmo em 2D)
    mesh->coords = malloc(mesh->n_nodes * 3 * sizeof(double));

	//Aloca a memoria do vetor de conectividades
    mesh->connectivity = malloc(mesh->n_elements * mesh->nodes_per_element * sizeof(int));

	//Verifica se a alocacao foi bem sucedida
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

    mesh_free_groups(mesh);
}

void mesh_free_groups(Mesh* mesh)
{
    for (int i = 0; i < mesh->n_groups; i++)
        free(mesh->groups[i].node_ids);

    free(mesh->groups);
    mesh->groups = NULL;
    mesh->n_groups = 0;
}

//Função que procura o grupo de nós pelo nome e retorna um ponteiro para ele, ou NULL se não encontrado. É usada para aplicar condições de contorno por nome, em vez de selecionar nós por coordenada.
//MeshGroup é um tipo criado na struct MeshGroup, que contém o nome do grupo, o número de nós e um ponteiro para os IDs dos nós
//groups[] é lido e criado na função read_groups_block() do arquivo unv_reader.c, durante a leitura do arquivo UNV
MeshGroup* mesh_find_group(Mesh* mesh, const char* name)
{
    for (int i = 0; i < mesh->n_groups; i++)
        if (strcmp(mesh->groups[i].name, name) == 0)
            return &mesh->groups[i];

    return NULL;
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

    printf("Groups: %d\n", mesh->n_groups);
    for (int i = 0; i < mesh->n_groups; i++)
        printf("  - %s (%d nos)\n", mesh->groups[i].name, mesh->groups[i].n_nodes);
}