#include "elasticity_2D.h"

void element_routine_2d(
    double* Ke,
    double* Re,
    Mesh* mesh,
    ElementType* etype,
    Material* mat,
    int e,
    double* u)
{
    double* gp_w = malloc(etype->n_gauss * 3 * sizeof(double));
    double* N = malloc(etype->nodes * sizeof(double));
    double* dN = malloc(etype->nodes * 2 * sizeof(double));

    etype->gauss_points(gp_w);

    for (int gp = 0; gp < etype->n_gauss; gp++)
    {
        /*
        Eu não preciso fazer:
        double* xi_ptr = &gp_w[gp * 3 + 0];
		double xi = *xi_ptr;

		Do jeito proposto, já é equivalente a acessar o valor de xi diretamente, sem a necessidade de criar um ponteiro intermediário.
        */
        double xi = gp_w[gp * 3 + 0];
        double eta = gp_w[gp * 3 + 1];
        double w = gp_w[gp * 3 + 2];

        /* 1 shape functions */
        etype->shape_functions(N, xi, eta);

        /* 2 derivatives */
        etype->shape_derivatives(dN, xi, eta);

        /* 3 compute Jacobian */
        compute_jacobian();

        /* 4 compute B matrix */
        compute_B();

        /* 5 strain */
        strain = B * u;

        /* 6 material model */
        mat->stress(sigma, C, strain);

        /* 7 stiffness */
        Ke += B ^ T * C * B * detJ * w;

        /* 8 residual */
        Re += B ^ T * sigma * detJ * w;
    }

    /* depois do loop de pontos de Gauss */
    free(N);
    free(dN);
    free(gp_w);
}