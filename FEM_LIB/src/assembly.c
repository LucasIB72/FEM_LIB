#include "assembly.h"
#include "element_routine.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Lista dinamica de pares (coluna, valor) para uma linha do CSR. A busca
   linear em row_entries_add() garante que cada coluna apareca no maximo
   uma vez por linha, independente da ordem em que os elementos que
   compartilham essa linha sao percorridos. */
typedef struct { int* col; double* val; int size, cap; } RowEntries;

static void row_entries_add(RowEntries* row, int col, double val)
{
    for (int k = 0; k < row->size; k++)
    {
        if (row->col[k] == col)
        {
            row->val[k] += val;
            return;
        }
    }

    if (row->size == row->cap)
    {
        row->cap = row->cap ? row->cap * 2 : 8;
        row->col = realloc(row->col, row->cap * sizeof(int));
        row->val = realloc(row->val, row->cap * sizeof(double));
    }

    row->col[row->size] = col;
    row->val[row->size] = val;
    row->size++;
}

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
       Monta, para cada linha (grau de liberdade), a lista de pares
       (coluna, valor) ja deduplicados e acumulados via row_entries_add().
       --------------------------------------------------------------- */

    RowEntries* rows = calloc(n, sizeof(RowEntries));

	//Percorre todos os elementos da malha
    for (int e = 0; e < mesh->n_elements; e++)
    {
        integrate_stiffness_matrix(Ke, Re, etype, physics, mesh, e, material_properties, NULL);

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
                        double val = Ke[(a * dof_per_node + i) * ndof_e +
                                        (b * dof_per_node + j)];

                        row_entries_add(&rows[gi], gj, val);
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
        RI[i + 1] = RI[i] + rows[i].size;

	//O ultimo elemento de RI indica o total de elementos não nulos na matriz de rigidez global. Ele seria o ultimo elemento do vetor "pointerE" da matriz esparsa no formato CSR. Então na verdade RI[0:n-1] = pointerB[0:n-1] e RI[n] = pointerE[n-1] em que n é o número de´graus de liberdade. Se usa n-1 nos vetores pois em C os indices começam em 0.
    int total_nnz = RI[n];

    int* COL = malloc(total_nnz * sizeof(int));
    double* VAL = malloc(total_nnz * sizeof(double));

    /* Concatena as listas por linha nos vetores finais COL/VAL do CSR */
    for (int i = 0; i < n; i++)
    {
        memcpy(COL + RI[i], rows[i].col, rows[i].size * sizeof(int));
        memcpy(VAL + RI[i], rows[i].val, rows[i].size * sizeof(double));
        free(rows[i].col);
        free(rows[i].val);
    }
    free(rows);

    /* ---------------------------------------------------------------
       Ordena as colunas de cada linha em ordem crescente. A insercao
       acima segue a ordem de varredura dos elementos/nos locais, nao a
       ordem das colunas globais - solvers diretos como o PARDISO exigem
       indices de coluna ordenados por linha no formato CSR
       --------------------------------------------------------------- */
    for (int i = 0; i < n; i++)
    {
        for (int p = RI[i] + 1; p < RI[i + 1]; p++)
        {
            int col_key = COL[p];
            double val_key = VAL[p];
            int q = p - 1;
            while (q >= RI[i] && COL[q] > col_key)
            {
                COL[q + 1] = COL[q];
                VAL[q + 1] = VAL[q];
                q--;
            }
            COL[q + 1] = col_key;
            VAL[q + 1] = val_key;
        }
    }

    /* ---------------------------------------------------------------
       Sa�das
       --------------------------------------------------------------- */

    *rowIndex = RI;
    *columns  = COL;
    *values   = VAL;
    *nnz      = total_nnz;

    free(Re);
    free(Ke);
}
