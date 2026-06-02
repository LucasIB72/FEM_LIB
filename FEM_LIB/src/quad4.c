#include "quad4.h"

void quad4_shape_functions(double* N, const double* natural_coords)
{
    N[0] = 0.25 * (1.0 - natural_coords[0]) * (1.0 - natural_coords[1]);
    N[1] = 0.25 * (1.0 + natural_coords[2]) * (1.0 - natural_coords[3]);
    N[2] = 0.25 * (1.0 + natural_coords[4]) * (1.0 + natural_coords[5]);
    N[3] = 0.25 * (1.0 - natural_coords[6]) * (1.0 + natural_coords[7]);
}


// dN: derivadas das funções de forma, passadas como referência para a função (o array deve ser alocado antes de chamar a função)
// natural_coords: coordenadas naturais (xi, eta, chi) do ponto onde as derivadas estão sendo calculadas, passadas como referência para a função, apenas como leitura (const)
void quad4_shape_derivatives(double* dN, const double* natural_coords)
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

ElementType create_quad4_element()
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