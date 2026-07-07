#include "bc.h"
#include <stdlib.h>
#include <string.h>

void apply_displacement_bc_dense(double* K, double* F, int n, const DisplacementBC* bc)
{
	//Loop sobre o numero de condições de contorno de deslocamento
    for (int b = 0; b < bc->n_dofs; b++)
    {
		//Obtem o grau de liberdade "d" e o valor do deslocamento "u_p" da condição de contorno
        int d = bc->dof[b];
        double u_p = bc->value[b];

        for (int j = 0; j < n; j++)
        {
            //Percorre linha a linha da matriz de rigidez com "n" graus de liberdade, e pega o valor K[j * n + d] da coluna "d" respectivo ao grau de liberdade prescrito
            if (j == d) continue;
            F[j] -= K[j * n + d] * u_p;
            K[j * n + d] = 0.0;
            K[d * n + j] = 0.0;
        }
		//Define o valor da diagonal K_dd como 1.0, e o valor do vetor de forças F[d] como o valor do deslocamento prescrito u_p
        K[d * n + d] = 1.0;
        F[d] = u_p;
    }
}

void apply_displacement_bc_csr(int* rowIndex, int* columns, double* values,
    double* F, int n, const DisplacementBC* bc)
{
    int* diag_pos = malloc(n * sizeof(int));
    if (!diag_pos) return;
    for (int i = 0; i < n; i++) diag_pos[i] = -1;

    for (int i = 0; i < n; i++)
        for (int p = rowIndex[i]; p < rowIndex[i + 1]; p++)
            if (columns[p] == i)
                { diag_pos[i] = p; break; }

    for (int b = 0; b < bc->n_dofs; b++)
    {
        int d = bc->dof[b];
        double u_p = bc->value[b];

        for (int p = rowIndex[d]; p < rowIndex[d + 1]; p++)
        {
            int j = columns[p];
            if (j == d) continue;
            F[j] -= values[p] * u_p;
        }

        for (int i = 0; i < n; i++)
        {
            if (i == d) continue;
            for (int p = rowIndex[i]; p < rowIndex[i + 1]; p++)
            {
                if (columns[p] == d)
                {
                    F[i] -= values[p] * u_p;
                    values[p] = 0.0;
                    break;
                }
            }
        }

        for (int p = rowIndex[d]; p < rowIndex[d + 1]; p++)
        {
            if (columns[p] != d)
                values[p] = 0.0;
        }

        if (diag_pos[d] >= 0)
            values[diag_pos[d]] = 1.0;

        F[d] = u_p;
    }

    free(diag_pos);
}

void apply_force_bc(double* F, const ForceBC* bc)
{
    for (int b = 0; b < bc->n_dofs; b++)
        F[bc->dof[b]] += bc->value[b];
}
