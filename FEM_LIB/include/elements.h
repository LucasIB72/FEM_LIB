#ifndef ELEMENTS_H
#define ELEMENTS_H

/**
* @file elements.h
* @brief Define a estrutura ElementType, que representa um tipo de elemento finito.
* @details A estrutura ElementType contém informações sobre o tipo de elemento, como a dimensão, o número de nós, os graus de liberdade por nó, o número de pontos de Gauss e ponteiros para funções que calculam as funções de forma, suas derivadas e os pontos de Gauss. ElementType contém o que é genérico de todos os tipos de elemento, e vai instanciar elementos específicos (ex: quad4). Portanto o sourcecode .c de cada tipo de elemento (ex: quad4.c) vai possuir uma função que retorna um ElementType preenchido com as informações específicas daquele tipo de elemento.
*/

typedef struct ElementType
{
    int dim;
    int nodes;
    int dof_per_node;
    int n_gauss;
    //Apontam para funcoes de forma de um elemento específico (ex: QUAD4)
    void (*shape_functions)(double* N, const double* xi);/**< É um ponteiro para uma função do elemento que calcula as funções de forma (shape functions) para um elemento finito. */
    void (*shape_derivatives)(double* dN, const double* xi);/**< É um ponteiro para uma função do elemento que calcula as derivadas das funções de forma para um elemento finito. */
	void (*gauss_points)(double* gp_w);/**< É um ponteiro para uma função do elemento que calcula os pontos de Gauss e seus pesos para um elemento finito. */
	//Aponta para uma física específica (ex: Elasticidade 2D)
    void (*integrate_physics)(double* Ke, double* Re, const double* dN_global, double detJ, double w, const double* material_properties, int n_nodes, int n_dim, int ndof, const double* u_e); /**< É um ponteiro para uma função da física escolhida, que calcula a parte física da integração das matrizes do elemento */
    

} ElementType;

#endif
