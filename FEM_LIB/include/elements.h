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
    const double* gauss_points;/**< É um ponteiro que aponta para um array que contém os pesos para cada ponto de Gauss do elemento. */
    
    void (*shape_functions)(double* N, const double* xi);/**< É um ponteiro para uma função que calcula as funções de forma (shape functions) para um elemento finito. */
    void (*shape_derivatives)(double* dN, const double* xi);/**< É um ponteiro para um array onde as funções de forma calculadas serão armazenadas. */
	void (*gauss_points)(double* gp_w);/**< É um ponteiro para uma função que calcula os pontos de Gauss e seus pesos para um elemento finito. */
    

} ElementType;

#endif
