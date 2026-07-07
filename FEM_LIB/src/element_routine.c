#include "element_routine.h"
#include "elements.h"
#include <stdlib.h>
#include <math.h>

// A funcao calcula a matriz de rigidez do elemento (Ke). Primeiro extrai as coordenadas nodais e, 
// se houver, os deslocamentos do elemento a partir do mesh. Em seguida processa a parte geometrica 
// da integracao (derivadas, Jacobiano) e delega a parte fisica para physics->integrate.
void integrate_stiffness_matrix(double* Ke, double* Re, ElementType* element, PhysicsModel* physics, Mesh* mesh, int e, const double* material_properties, const double* u)
{
    int n_nodes = element->nodes;
    int n_gp = element->n_gauss;
	int n_dim = element->dim;
    int ndof = n_nodes * physics->dof_per_node;

    // Extrai as coordenadas nodais do elemento "e" a partir do mesh
    double* node_coords = malloc(n_nodes * 3 * sizeof(double));
    if (!node_coords) return;

    for (int j = 0; j < n_nodes; j++)
    {
        int node_id = mesh->connectivity[e * mesh->nodes_per_element + j];
        node_coords[j * 3 + 0] = mesh->coords[node_id * 3 + 0];
        node_coords[j * 3 + 1] = mesh->coords[node_id * 3 + 1];
        node_coords[j * 3 + 2] = mesh->coords[node_id * 3 + 2];
    }

	// Extrai o vetor de deslocamentos do elemento "e" a partir do vetor global "u" (se tiver sido calculado antes)
    double* u_e = NULL;
    if (u)
    {
        u_e = malloc(ndof * sizeof(double));
        if (!u_e) { free(node_coords); return; }
    }
    if (u)
    {
        for (int j = 0; j < n_nodes; j++)
        {
            int node_id = mesh->connectivity[e * mesh->nodes_per_element + j];
            for (int d = 0; d < physics->dof_per_node; d++)
            {
                u_e[j * physics->dof_per_node + d] = u[node_id * physics->dof_per_node + d];
            }
        }
    }
        

	// Aloca memoria para as derivadas das funcoes de forma e os pontos de Gauss com pesos
    double* dN = malloc(n_dim * n_nodes * sizeof(double));
    double* gp_w = malloc((n_dim + 1) * n_gp * sizeof(double));
	double* dN_global = malloc(n_dim * n_nodes * sizeof(double));
    double* xi_arr = malloc(n_dim * sizeof(double));

    // Procura erro de alocacao
    if (!dN || !gp_w || !dN_global || !xi_arr)
    {
        free(u_e);
        free(node_coords);
        free(dN);
        free(gp_w);
        free(dN_global);
        free(xi_arr);
        return;
    }

	// Chama a funcao do elemento para obter os pontos de Gauss e seus pesos
    element->gauss_points(gp_w);

    // Zera as matrizes do elemento
    for (int i = 0; i < ndof * ndof; i++)
        Ke[i] = 0.0;
    for (int i = 0; i < ndof; i++)
        Re[i] = 0.0;

	// Loop sobre os pontos de Gauss
    for (int i = 0; i < n_gp; i++)
    {
		// Extrai as coordenadas naturais e o peso do ponto de Gauss 
		// xi_arr = {xi_1, eta_1, zeta_1, xi_2, eta_2, zeta_2, ...} e w = {w_1, w_2, ...}
        for (int d = 0; d < n_dim; d++)
            xi_arr[d] = gp_w[i * (n_dim + 1) + d];
        double w = gp_w[i * (n_dim + 1) + n_dim];

		// Calcula as derivadas das funcoes de forma no ponto de Gauss, em relacao as coordenadas naturais
        element->shape_derivatives(dN, xi_arr);

		// Calcula a matriz Jacobiana J = {{dx/dxi, dx/deta, dx/dzeta}, {dy/dxi, dy/deta, dy/dzeta}, {dz/dxi, dz/deta, dz/dzeta}}
        double J[9] = {0.0};
        for (int j = 0; j < n_dim; j++)
            for (int i = 0; i < n_dim; i++)
                for (int k = 0; k < n_nodes; k++)
                    J[i * n_dim + j] += dN[k * n_dim + i] * node_coords[k * 3 + j];

        // Calcula o determinante da matriz Jacobiana e sua inversa J_inv.
        double detJ;
        double J_inv[9];
        if (n_dim == 2)
        {
            detJ = J[0] * J[3] - J[1] * J[2];
            J_inv[0] =  J[3] / detJ;
            J_inv[1] = -J[1] / detJ;
            J_inv[2] = -J[2] / detJ;
            J_inv[3] =  J[0] / detJ;
        }
        else
        {
            detJ = J[0] * (J[4] * J[8] - J[5] * J[7])
                 - J[1] * (J[3] * J[8] - J[5] * J[6])
                 + J[2] * (J[3] * J[7] - J[4] * J[6]);
            double inv_det = 1.0 / detJ;
            J_inv[0] = (J[4] * J[8] - J[5] * J[7]) * inv_det;
            J_inv[1] = (J[2] * J[7] - J[1] * J[8]) * inv_det;
            J_inv[2] = (J[1] * J[5] - J[2] * J[4]) * inv_det;
            J_inv[3] = (J[5] * J[6] - J[3] * J[8]) * inv_det;
            J_inv[4] = (J[0] * J[8] - J[2] * J[6]) * inv_det;
            J_inv[5] = (J[2] * J[3] - J[0] * J[5]) * inv_det;
            J_inv[6] = (J[3] * J[7] - J[4] * J[6]) * inv_det;
            J_inv[7] = (J[1] * J[6] - J[0] * J[7]) * inv_det;
            J_inv[8] = (J[0] * J[4] - J[1] * J[3]) * inv_det;
        }

		// Calcula as derivadas das funcoes de forma em relacao as coordenadas globais dN_global = {dN1/dx, dN1/dy, dN2/dx, dN2/dy, ..., dNn/dx, dNn/dy}
		// dN = {dN1/dxi, dN1/deta, dN2/dxi, dN2/deta, ..., dNn/dxi, dNn/deta}
		// {dN1/dx, dN1/dy} = {dN1/dxi, dN1/deta} [J_inv[0] J_inv[1]
		//  							           J_inv[2] J_inv[3]]
        for (int node = 0; node < n_nodes; node++)
        {
            for (int i_dir_glob = 0; i_dir_glob < n_dim; i_dir_glob++)
            {
                dN_global[node * n_dim + i_dir_glob] = 0.0;
                for (int i_dir_loc = 0; i_dir_loc < n_dim; i_dir_loc++)
                    dN_global[node * n_dim + i_dir_glob] += J_inv[i_dir_glob * n_dim + i_dir_loc] * dN[node * n_dim + i_dir_loc];
            }
		}

        // Chama a funcao de integracao fisica especifica do elemento
        if (physics && physics->integrate)
            physics->integrate(Ke, Re, dN_global, detJ, w, material_properties, n_nodes, n_dim, ndof, u_e);
    }

    free(u_e);
    free(node_coords);
    free(xi_arr);
    free(dN);
    free(gp_w);
    free(dN_global);
}