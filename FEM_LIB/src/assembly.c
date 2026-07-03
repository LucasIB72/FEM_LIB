#include "assembly.h"
#include <stdlib.h>
#include <string.h>

void assemble_global_stiffness(
    double* K_global,
    double* R_global,
    Mesh* mesh,
    ElementType* etype,
    Material* mat,
    const double* u)
{
    int n_nodes_e = etype->nodes;
    int dof_per_node = etype->dof_per_node;
    int ndof_e = n_nodes_e * dof_per_node;
    int total_dofs = mesh->total_dofs;

    double* Ke = malloc(ndof_e * ndof_e * sizeof(double));
    double* Re = malloc(ndof_e * sizeof(double));

    //Se nao conseguir alocar, sai da rotina
    if (!Ke || !Re)
    {
        free(Ke);
        free(Re);
        return;
    }


    for (int e = 0; e < mesh->n_elements; e++)
    {
        element_routine_e2d(Ke, Re, mesh, etype, mat, e, u);

        for (int a = 0; a < n_nodes_e; a++)
        {
            int ia = mesh->connectivity[e * mesh->nodes_per_element + a];
            for (int i = 0; i < dof_per_node; i++)
            {
                int gi = ia * dof_per_node + i;
                for (int b = 0; b < n_nodes_e; b++)
                {
                    int ib = mesh->connectivity[e * mesh->nodes_per_element + b];
                    for (int j = 0; j < dof_per_node; j++)
                    {
                        int gj = ib * dof_per_node + j;
                        K_global[gi + gj * total_dofs] +=
                            Ke[(a * dof_per_node + i) * ndof_e +
                               (b * dof_per_node + j)];
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
    Material* mat)
{
    int n_nodes_e = etype->nodes;
    int dof_per_node = etype->dof_per_node;
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
       PASS 1 : conta não-nulos únicos por linha.
       Usa o próprio gi como carimbo no marker[gj]: se marker[gj] == gi
       a coluna já foi contada para esta linha; caso contrário é nova.
       --------------------------------------------------------------- */

    int* nnz_per_row = calloc(n, sizeof(int));
    int* marker = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) marker[i] = -1;

    for (int e = 0; e < mesh->n_elements; e++)
    {
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
                        if (marker[gj] != gi)
                        {
                            marker[gj] = gi;
                            nnz_per_row[gi]++;
                        }
                    }
                }
            }
        }
    }

    /* ---------------------------------------------------------------
       Constrói rowIndex a partir da contagem
       --------------------------------------------------------------- */

    int* RI = malloc((n + 1) * sizeof(int));
    RI[0] = 0;
    for (int i = 0; i < n; i++)
        RI[i + 1] = RI[i] + nnz_per_row[i];

    int total_nnz = RI[n];

    int* COL = malloc(total_nnz * sizeof(int));
    double* VAL = calloc(total_nnz, sizeof(double));

    /* ---------------------------------------------------------------
       PASS 2 : preenche colunas e acumula valores.
       marker[gj] == gi  →  coluna já existe, busca e acumula
       marker[gj] != gi  →  coluna inédita (ou de outra linha), armazena
       --------------------------------------------------------------- */

    int* row_pos = malloc(n * sizeof(int));
    memcpy(row_pos, RI, n * sizeof(int));
    for (int i = 0; i < n; i++) marker[i] = -1;

    for (int e = 0; e < mesh->n_elements; e++)
    {
        element_routine_e2d(Ke, Re, mesh, etype, mat, e, NULL);

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
       Saídas
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
