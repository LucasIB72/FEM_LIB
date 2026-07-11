#ifndef MESH_H
#define MESH_H

/* Um grupo nomeado de nos da malha (ex: "engaste", "carga_topo"), tipicamente
   definido no pre-processador (SALOME) e usado para aplicar condicoes de
   contorno por nome, em vez de selecionar nos por coordenada. */
typedef struct MeshGroup
{
    char name[64];
    int n_nodes;
    int* node_ids;      // indices 0-based em mesh->coords
} MeshGroup;

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

    int n_groups;
    MeshGroup* groups;      // [n_groups]

} Mesh;

void mesh_init(Mesh* mesh);

void mesh_allocate(Mesh* mesh);

void mesh_free(Mesh* mesh);

void mesh_free_groups(Mesh* mesh);

MeshGroup* mesh_find_group(Mesh* mesh, const char* name);

void mesh_print_info(const Mesh* mesh);

#endif