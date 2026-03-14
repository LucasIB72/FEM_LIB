#include "quad4.h"

void quad4_shape_functions(double xi, double eta, double* N)
{
    N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
    N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
    N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
    N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
}

void quad4_shape_derivatives(double xi, double eta, double dN[4][2])
{
    dN[0][0] = -0.25 * (1.0 - eta);
    dN[0][1] = -0.25 * (1.0 - xi);

    dN[1][0] = 0.25 * (1.0 - eta);
    dN[1][1] = -0.25 * (1.0 + xi);

    dN[2][0] = 0.25 * (1.0 + eta);
    dN[2][1] = 0.25 * (1.0 + xi);

    dN[3][0] = -0.25 * (1.0 + eta);
    dN[3][1] = 0.25 * (1.0 - xi);
}