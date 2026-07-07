#include "elasticity_2D_plane_stress.h"
#include <stdlib.h>

void elasticity_2d_plane_stress_integrate(
    double* Ke, double* Re,
    const double* dN_global,
    double detJ, double w,
    const double* material_properties,
    int n_nodes, int n_dim, int ndof,
    const double* u_e)
{
    double E = material_properties[0];
    double nu = material_properties[1];

	// Matriz de elasticidade para o estado plano de tensoes
    double c = E / (1.0 - nu * nu);
    double C[3][3] = {
        {c, c * nu, 0.0},
        {c * nu, c, 0.0},
        {0.0, 0.0, c * (1.0 - nu) / 2.0}
    };

    double* B_data = calloc(3 * ndof, sizeof(double));
    if (!B_data) return;
    double* B[3] = {
        B_data + 0 * ndof,
        B_data + 1 * ndof,
        B_data + 2 * ndof
    };

    for (int j = 0; j < n_nodes; j++)
    {
        B[0][j * 2]     = dN_global[j * n_dim + 0];
        B[1][j * 2 + 1] = dN_global[j * n_dim + 1];
        B[2][j * 2]     = dN_global[j * n_dim + 1];
        B[2][j * 2 + 1] = dN_global[j * n_dim + 0];
    }

    double** K = malloc(ndof * sizeof(double*));
    if (!K) { free(B_data); return; }
    for (int i = 0; i < ndof; i++)
        K[i] = Ke + i * ndof;


	//Forma indicial do produto matricial [Be]^T * [C] * [Be]
    for (int a = 0; a < ndof; a++)
    {
        for (int b = 0; b < ndof; b++)
        {
            double sum = 0.0;
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                    sum += B[k][a] * C[k][l] * B[l][b];
            K[a][b] += sum * detJ * w;
        }
    }

	//Calcula as metricas pos processadas: deformacao, tensao e reacoes internas (residuos)
    if (u_e)
    {
        double strain[3] = { 0.0, 0.0, 0.0 };
        for (int k = 0; k < 3; k++)
            for (int a = 0; a < ndof; a++)
                strain[k] += B[k][a] * u_e[a];

        double sigma[3] = { 0.0, 0.0, 0.0 };
        for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++)
                sigma[k] += C[k][l] * strain[l];

        for (int a = 0; a < ndof; a++)
        {
            double sum = 0.0;
            for (int k = 0; k < 3; k++)
                sum += B[k][a] * sigma[k];
            Re[a] += sum * detJ * w;
        }
    }

    free(K);
    free(B_data);
}
