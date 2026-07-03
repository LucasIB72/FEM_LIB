#include "simp_input.h"
#include <stdlib.h>
#include <math.h>

int generate_rect_mesh(
    double L_x, double L_y, double e_s,
    double** coords,
    int** connectivity,
    int* n_nodes,
    int* n_elements,
    int* nodes_per_element)
{

	// Verifica se os parametros de dimensao do retangulo e o tamanho do elemento sao validos
    if (L_x <= 0.0 || L_y <= 0.0 || e_s <= 0.0)
        return -1;

	// Calcula o numero de elementos em cada direcao
    int nx = (int)(L_x / e_s);
    int ny = (int)(L_y / e_s);

	// Garante que haja pelo menos um elemento em cada direcao
    if (nx < 1) nx = 1;
    if (ny < 1) ny = 1;

    // Calcula o tamanho do elemento em cada direcao
    double dx = L_x / nx;
    double dy = L_y / ny;

	// Calcula as propriedades da malha e aloca nos enderecos de memoria apontado pelos ponteiros fornecidos
    *n_nodes = (nx + 1) * (ny + 1);
    *n_elements = nx * ny;
    *nodes_per_element = 4;

	// Aloca memoria para o vetor de coordenadas (x1, y1, z1, x2, y2, z2, ..., xn, yn, zn) para cada no da malha
    *coords = malloc(*n_nodes * 3 * sizeof(double));

	// Aloca memoria para o vetor de conectividades (apenas elementos quadrilaterais, cada elemento tem 4 nos) para cada elemento da malha
    *connectivity = malloc(*n_elements * 4 * sizeof(int));

    if (*coords == NULL || *connectivity == NULL)
    {
        free(*coords);
        free(*connectivity);
        return -1;
    }

	// Preenche o vetor de coordenadas a contagem comeca da esquerda para a direita e de baixo para cima, com z = 0.0
	int node_id = 0;
    for (int j = 0; j <= ny; j++)
    {
        for (int i = 0; i <= nx; i++)
        {
            
            (*coords)[node_id * 3 + 0] = i * dx;
            (*coords)[node_id * 3 + 1] = j * dy;
            (*coords)[node_id * 3 + 2] = 0.0;

            node_id++;
        }
    }

	// Preenche o vetor de conectividades considerando os nos no sentido anti-horario, comecando do no inferior esquerdo de cada elemento
	int elem_id = 0;
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            
            int row = j * (nx + 1);

            (*connectivity)[elem_id * 4 + 0] = row + i;
            (*connectivity)[elem_id * 4 + 1] = row + i + 1;
            (*connectivity)[elem_id * 4 + 2] = row + (nx + 1) + i + 1;
            (*connectivity)[elem_id * 4 + 3] = row + (nx + 1) + i;

            elem_id++;
        }
    }

    return 0;
}

int generate_rect_mesh_to_mesh(
    double L_x, double L_y, double e_s,
    Mesh* mesh)
{
    int ret = generate_rect_mesh(
        L_x, L_y, e_s,
        &mesh->coords,
        &mesh->connectivity,
        &mesh->n_nodes,
        &mesh->n_elements,
        &mesh->nodes_per_element);

    if (ret == 0)
    {
        mesh->dim = 2;
        mesh->dof_per_node = 2;
        mesh->total_dofs = mesh->n_nodes * mesh->dof_per_node;
    }

    return ret;
}
