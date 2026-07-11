#include "out_writter.h"
#include <stdio.h>
#include <math.h>

void write_vtk(const char* filename, Mesh* mesh, double* u)
{
    //O arquivo VTK (Visualization Toolkit) é o padrão para visualização de dados científicos no paraview.
	//O arquivo VTK suporta diferentes tipos de malhas, estruturadas ou não: uniform rectilinear grid, rectilinear grid, curvilinear grid, AMR dataset e unstructured grid.
    //O formato unstructured grid é o mais generalista, mas ocupa mais memória
    //POINTS: define os pontos da malha (coordenadas dos nós)
    //CELLS: são entidades geométricas que definem regiões no espaço (são os elementos)
	//CELL_TYPES: define o tipo de elemento (ex: VTK_QUAD, VTK_TRIANGLE, VTK_TETRA, etc). O tipo é indicado por um índice
	//POINT_DATA: define os dados associados aos pontos da malha (ex: deslocamentos, temperatura, etc)
	//CELL_DATA: define os dados associados aos elementos da malha (ex: tensões, deformações, etc)
    //
	// Dentro de POINT_DATA e CELL_DATA estão escritos os "attributes blocks", que podem ser campos escalares (SCALARS), vetoriais (VECTORS), tensoriais (TENSORS), campo genérico (FIELD), etc.
    //
	// Estrutura do arquivo: header (sempre 4 linhas) + data sections (indroduzidas por keywords)
    //
    // Header:
	// # vtk DataFile Version 3.0 (linha 1): fala o formado do arquivo vtk, que é a versão 3.0 legacy
	// FEM_LIB (linha 2): é livre para o usuário colocar qualquer comentário
	// ASCII (linha 3): fala o formato dos dados, que pode ser ASCII ou BINARY
	// DATASET UNSTRUCTURED_GRID (linha 4): define a topologia dos dados
    //
    // Data sections:
	// POINTS n double: define os pontos da malha, onde n é o número de nós e double é o tipo de dado (float ou double)
	// 0.000000e+00 0.000000e+00 0.000000e+00 (coordenadas x, y, z do nó 0)
	// 0.000000e+00 1.000000e+00 0.000000e+00 (coordenadas x, y, z do nó 1)
    // ...
	// 1.000000e+00 1.000000e+00 0.000000e+00 (coordenadas x, y, z do nó n-1)
    //
	// CELLS n m: define os elementos da malha, onde n é o número de elementos e m é o número total de índices (número de elementos * (número de nós por elemento + 1))
	// n_nodes_per_element node0 node1 node2 node3 (elem 1)
	// n_nodes_per_element node0 node1 node2 node3 (elem 2)
    // ...
	// n_nodes_per_element node0 node1 node2 node3 (elem n)
    //
	// CELL_TYPES n: define o tipo de elemento
	// VTK_ELEMENT_TYPE (elem 1)
	// VTK_ELEMENT_TYPE (elem 2)
    // ...
	// VTK_ELEMENT_TYPE (elem n)
    //
	// POINT_DATA n: define os dados associados aos pontos da malha, onde n é o número de nós
	// VECTORS field_name double: define um campo vetorial com o nome field_name e tipo double
	// ux uy uz (deslocamento no nó 0)
	// ux uy uz (deslocamento no nó 1)
    // ...
	// ux uy uz (deslocamento no nó n-1)
	// SCALARS field_name data_type ncomp: define um campo escalar com o nome field_name o tipo data_type e ncomp componentes (ncomp = 1 para campo escalar, 2 vetor 2 componentes, até 4)
	// LOOKUP_TABLE default: define a tabela de cores padrão para o plot
	// value (campo escalar no nó 0)
	// value (campo escalar no nó 1)
    // ...
	// value (campo escalar no nó n-1)

    FILE* fp = fopen(filename, "w");
    if (!fp) return;

    fprintf(fp, "# vtk DataFile Version 3.0\nFEM_LIB\nASCII\nDATASET UNSTRUCTURED_GRID\n");
    fprintf(fp, "POINTS %d double\n", mesh->n_nodes);
    for (int i = 0; i < mesh->n_nodes; i++)
        fprintf(fp, "%e %e %e\n", mesh->coords[i * 3], mesh->coords[i * 3 + 1], mesh->coords[i * 3 + 2]);

    fprintf(fp, "CELLS %d %d\n", mesh->n_elements, mesh->n_elements * (mesh->nodes_per_element + 1));
    for (int e = 0; e < mesh->n_elements; e++)
    {
        fprintf(fp, "%d", mesh->nodes_per_element);
        for (int j = 0; j < mesh->nodes_per_element; j++)
            fprintf(fp, " %d", mesh->connectivity[e * mesh->nodes_per_element + j]);
        fprintf(fp, "\n");
    }

    fprintf(fp, "CELL_TYPES %d\n", mesh->n_elements);
    for (int e = 0; e < mesh->n_elements; e++)
        fprintf(fp, "%d\n", mesh->nodes_per_element == 4 ? 9 : 5);

    if (u)
    {
        fprintf(fp, "POINT_DATA %d\n", mesh->n_nodes);

        /* Campo vetorial (usado no Warp By Vector) */
        fprintf(fp, "VECTORS displacement double\n");
        for (int i = 0; i < mesh->n_nodes; i++)
            fprintf(fp, "%e %e 0.0\n", u[i * 2], u[i * 2 + 1]);

        /* Campo escalar de magnitude (facilita a coloracao direta) */
        fprintf(fp, "SCALARS displacement_magnitude double 1\n");
        fprintf(fp, "LOOKUP_TABLE default\n");
        for (int i = 0; i < mesh->n_nodes; i++)
        {
            double ux = u[i * 2], uy = u[i * 2 + 1];
            fprintf(fp, "%e\n", sqrt(ux * ux + uy * uy));
        }
    }

    fclose(fp);
}
