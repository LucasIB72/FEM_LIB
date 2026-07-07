#include "quad4.h"

void quad4_shape_functions(double* N, const double* xi)
{
    double eta_val = xi[1];
    double xi_val = xi[0];

    N[0] = 0.25 * (1.0 - xi_val) * (1.0 - eta_val);
    N[1] = 0.25 * (1.0 + xi_val) * (1.0 - eta_val);
    N[2] = 0.25 * (1.0 + xi_val) * (1.0 + eta_val);
    N[3] = 0.25 * (1.0 - xi_val) * (1.0 + eta_val);
}


// A função quad4_shape_derivatives calcula as derivadas das funções de forma em relação às coordenadas naturais (xi e eta) para um elemento quadrilateral de 4 nós. As derivadas são armazenadas no array dN, onde as primeiras 4 posições correspondem às derivadas em relação a xi e as próximas 4 posições correspondem às derivadas em relação a eta.
// dN: derivadas das funções de forma, passadas como referência para a função (o array deve ser alocado antes de chamar a função)
// eta e xi: coordenadas naturais, passadas como referência para a função.
void quad4_shape_derivatives(double* dN, const double* xi)
{

	double eta_val = xi[1];
	double xi_val = xi[0];

    // Derivadas da função de forma N_1
	dN[0] = -0.25 * (1.0 - eta_val); //dN_1/dxi
	dN[1] = -0.25 * (1.0 - xi_val); //dN_1/deta

    // Derivadas da função de forma N_2
	dN[2] = 0.25 * (1.0 - eta_val); //dN_2/dxi
	dN[3] = -0.25 * (1.0 + xi_val); //dN_2/deta

    // Derivadas da função de forma N_3
	dN[4] = 0.25 * (1.0 + eta_val); //dN_3/dxi
	dN[5] = 0.25 * (1.0 + xi_val); //dN_3/deta

    // Derivadas da função de forma N_4
	dN[6] = -0.25 * (1.0 + eta_val); //dN_4/dxi
	dN[7] = 0.25 * (1.0 - xi_val); //dN_4/deta
}

void quad4_gauss_points(double* gp_w)
{
	gp_w[0] = -0.5773502691896257; gp_w[1] = -0.5773502691896257; gp_w[2] = 1.0;
	gp_w[3] = 0.5773502691896257; gp_w[4] = -0.5773502691896257; gp_w[5] = 1.0;        
	gp_w[6] = 0.5773502691896257; gp_w[7] = 0.5773502691896257; gp_w[8] = 1.0; 
	gp_w[9] = -0.5773502691896257; gp_w[10] = 0.5773502691896257; gp_w[11] = 1.0;    
}

ElementType create_quad4_element()
{
    ElementType e;

    e.dim = 2;
    e.nodes = 4;
    //e.dof_per_node = 2; // Movido para PhysicsModel
    e.n_gauss = 4;

    e.shape_functions = quad4_shape_functions;
    e.shape_derivatives = quad4_shape_derivatives;
	e.gauss_points = quad4_gauss_points;

    return e;
}