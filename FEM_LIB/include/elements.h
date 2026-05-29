#ifndef ELEMENTS_H
#define ELEMENTS_H

typedef struct ElementType
{
    int dim;
    int nodes;
    int dof_per_node;
    int n_gauss;
    
	//*shape_functions é um ponteiro para uma função que calcula as funções de forma (shape functions) para um elemento finito.
	//double* N é um ponteiro para um array onde as funções de forma calculadas serão armazenadas.
    void (*shape_functions)(double* N, double xi, double eta);
    void (*shape_derivatives)(double* dN, double xi, double eta);
    void (*gauss_points)(double* gp_w);

} ElementType;

#endif
