#include "assembly.h"
#include "element_routine.h"
#include <stdlib.h>
#include <string.h>

void assemble_global_stiffness(
    double* K_global,
    double* R_global,
    Mesh* mesh,
    ElementType* etype,
    PhysicsModel* physics,
    const double* material_properties,
    const double* u)
{
    int n_nodes_e = etype->nodes;
    int dof_per_node = physics->dof_per_node;
    int ndof_e = n_nodes_e * dof_per_node;
    int total_dofs = mesh->total_dofs;

    double* Ke = malloc(ndof_e * ndof_e * sizeof(double));
    double* Re = malloc(ndof_e * sizeof(double));

    if (!Ke || !Re)
    {
        free(Ke);
        free(Re);
        return;
    }

    for (int e = 0; e < mesh->n_elements; e++)
    {
		//Obtem a matriz de rigidez do elemento e o vetor de reacoes internas (residuos) do elemento
        integrate_stiffness_matrix(Ke, Re, etype, physics, mesh, e, material_properties, u);

        for (int a = 0; a < n_nodes_e; a++)
        {
			//Obtem o indice global "ia" correspondente ao no "a" do elemento "e"
            int ia = mesh->connectivity[e * mesh->nodes_per_element + a];
            for (int i = 0; i < dof_per_node; i++)
            {
				//Calcula o indice global "gi" correspondente ao grau de liberdade "i" do no "a" do elemento "e"
                int gi = ia * dof_per_node + i;
				//Calcula o indice elementar "gi_e" correspondente ao grau de liberdade "i" do no "a" do elemento "e"
				int gi_e = a * dof_per_node + i;
                for (int b = 0; b < n_nodes_e; b++)
                {
					//Obtem o indice do no global "ib" correspondente ao no "b" do elemento "e"
                    int ib = mesh->connectivity[e * mesh->nodes_per_element + b];
                    for (int j = 0; j < dof_per_node; j++)
                    {
						//Obtem o indice global "gj" correspondente ao grau de liberdade "j" do no "b" do elemento "e"
                        int gj = ib * dof_per_node + j;
						//Calcula o indice elementar "gj_e" correspondente ao grau de liberdade "j" do no "b" do elemento "e"
						int gj_e = b * dof_per_node + j;

						//Acumula o valor da matriz de rigidez do elemento na matriz de rigidez global
                        K_global[gi * total_dofs + gj] += Ke[gi_e * ndof_e + gj_e];
                    }
                }
                if (u)
                    R_global[gi] += Re[a * dof_per_node + i];
            }
        }
    }

    free(Re);
    free(Ke);
}

void assemble_global_stiffness_sparse(
    int** rowIndex,
    int** columns,
    double** values,
    int* nnz,
    Mesh* mesh,
    ElementType* etype,
    PhysicsModel* physics,
    const double* material_properties)
{
    int n_nodes_e = etype->nodes;
    int dof_per_node = physics->dof_per_node;
    int ndof_e = n_nodes_e * dof_per_node;
    int n = mesh->total_dofs;

    double* Ke = malloc(ndof_e * ndof_e * sizeof(double));
    double* Re = malloc(ndof_e * sizeof(double));

    if (!Ke || !Re)
    {
        free(Ke);
        free(Re);
        return;
    }

    /* ---------------------------------------------------------------
       PASS 1 : conta nao nulos unicos por linha.
       Usa o proprio gi como carimbo no marker[gj]: se marker[gj] == gi
       a coluna j� foi contada para esta linha; caso contr�rio � nova.
       --------------------------------------------------------------- */

    int* nnz_per_row = calloc(n, sizeof(int));
    int* marker = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) marker[i] = -1;

	//Percorre todos os elementos da malha
    for (int e = 0; e < mesh->n_elements; e++)
    {
        //Obtem a conectivade do elemento "e"
        int* conn_e = &mesh->connectivity[e * mesh->nodes_per_element];

		//Percorre todos os nos "a" do elemento "e"
        for (int a = 0; a < n_nodes_e; a++)
        {
			//Percorre os graus de liberdade de cada nó "a" do elemento "e"
            for (int i = 0; i < dof_per_node; i++)
            {
				//Calcula o indice global "gi" correspondente ao grau de liberdade "i" do no "a" do elemento "e"
                int gi = conn_e[a] * dof_per_node + i;

				//Percorre todos os nos "b" do elemento "e"
                for (int b = 0; b < n_nodes_e; b++)
                {
					//Percorre os graus de liberdade de cada nó "b" do elemento "e"
                    for (int j = 0; j < dof_per_node; j++)
                    {
						//Calcula o indice global "gj" correspondente ao grau de liberdade "j" do no "b" do elemento "e"
                        int gj = conn_e[b] * dof_per_node + j;
                        if (marker[gj] != gi)
                        {
                            marker[gj] = gi;
							//Calcula o numero nao nulos para a linha "gi" da matriz de rigidez global
                            nnz_per_row[gi]++;
                        }
                    }
                }
            }
        }
    }

   
	//O formato CSR pode trabalhar com os vetores que indical o inicio e o fim de cada linha da matriz esparsa. O vetor "pointerB" indica o índice do inicio de cada linha da matriz K no vetor de valores, e o vetor "pointerE" indica índice do fim de cada linha. Ou pode trabalhar com um unico vetor chamado "rowIndex" que indica o inicio de cada linha, e o fim de cada linha é dado pelo inicio da proxima linha. O vetor "rowIndex" tem tamanho n+1, onde n é o numero de graus de liberdade. O ultimo elemento do vetor "rowIndex" indica o total de elementos nao nulos na matriz esparsa.
    
	//RI é o vetor "rowIndex"
    int* RI = malloc((n + 1) * sizeof(int));
    RI[0] = 0;
    for (int i = 0; i < n; i++)
        RI[i + 1] = RI[i] + nnz_per_row[i];

	//O ultimo elemento de RI indica o total de elementos não nulos na matriz de rigidez global. Ele seria o ultimo elemento do vetor "pointerE" da matriz esparsa no formato CSR. Então na verdade RI[0:n-1] = pointerB[0:n-1] e RI[n] = pointerE[n-1] em que n é o número de´graus de liberdade. Se usa n-1 nos vetores pois em C os indices começam em 0.
    int total_nnz = RI[n];

    int* COL = malloc(total_nnz * sizeof(int));
    double* VAL = calloc(total_nnz, sizeof(double));

	//Define row_pos como um vetor que vai armazenar a posição atual de cada linha na matriz esparsa. Ele é inicializado com os valores de RI, que indicam o início de cada linha. À medida que os elementos não nulos são adicionados à matriz esparsa, row_pos é incrementado para apontar para a próxima posição disponível na linha correspondente.
    int* row_pos = malloc(n * sizeof(int));
    memcpy(row_pos, RI, n * sizeof(int));

    for (int i = 0; i < n; i++) marker[i] = -1;

    for (int e = 0; e < mesh->n_elements; e++)
    {
        integrate_stiffness_matrix(Ke, Re, etype, physics, mesh, e, material_properties, NULL);

        int* conn_e = &mesh->connectivity[e * mesh->nodes_per_element];

        for (int a = 0; a < n_nodes_e; a++)
        {
            for (int i = 0; i < dof_per_node; i++)
            {
                int gi = conn_e[a] * dof_per_node + i;

                for (int b = 0; b < n_nodes_e; b++)
                {
                    for (int j = 0; j < dof_per_node; j++)
                    {
                        int gj = conn_e[b] * dof_per_node + j;
                        double val = Ke[(a * dof_per_node + i) * ndof_e +
                                        (b * dof_per_node + j)];

                        if (marker[gj] != gi)
                        {
                            marker[gj] = gi;
                            COL[row_pos[gi]] = gj;
                            VAL[row_pos[gi]] = val;
                            row_pos[gi]++;
                        }
                        else
                        {
                            for (int p = RI[gi]; p < row_pos[gi]; p++)
                            {
                                if (COL[p] == gj)
                                {
                                    VAL[p] += val;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /* ---------------------------------------------------------------
       Sa�das
       --------------------------------------------------------------- */

    *rowIndex = RI;
    *columns  = COL;
    *values   = VAL;
    *nnz      = total_nnz;

    free(row_pos);
    free(marker);
    free(nnz_per_row);
    free(Re);
    free(Ke);
}
