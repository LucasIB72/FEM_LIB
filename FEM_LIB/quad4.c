#include "quad4.h"

void quad4_shape_functions(double* N, double xi, double eta)
{
    N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
    N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
    N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
    N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
}

void quad4_shape_derivatives(double* dN, double xi, double eta)
{
    dN[0] = -0.25 * (1.0 - eta);
    dN[4] = -0.25 * (1.0 - xi);

    dN[1] = 0.25 * (1.0 - eta);
    dN[5] = -0.25 * (1.0 + xi);

    dN[2] = 0.25 * (1.0 + eta);
    dN[6] = 0.25 * (1.0 + xi);

    dN[3] = -0.25 * (1.0 + eta);
    dN[7] = 0.25 * (1.0 - xi);
}

void quad4_gauss_points(double* gp_w)
{
	gp_w[0] = -0.5773502691896257; gp_w[4] = -0.5773502691896257; gp_w[8] = 1.0;
	gp_w[1] = 0.5773502691896257; gp_w[5] = -0.5773502691896257; gp_w[9] = 1.0;        
	gp_w[2] = 0.5773502691896257; gp_w[6] = 0.5773502691896257; gp_w[10] = 1.0; 
	gp_w[3] = -0.5773502691896257; gp_w[7] = 0.5773502691896257; gp_w[11] = 1.0;    
}

ElementType create_q4_element()
{
    ElementType e;

    e.dim = 2;
    e.nodes = 4;
    e.dof_per_node = 2;
    e.n_gauss = 4;

    e.shape_functions = quad4_shape_functions;
    e.shape_derivatives = quad4_shape_derivatives;
	e.gauss_points = quad4_gauss_points;

    return e;
}