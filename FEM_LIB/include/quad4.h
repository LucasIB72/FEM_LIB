#ifndef QUAD4_H
#define QUAD4_H

#include "elements.h"

void quad4_shape_functions(double* N, const double* natural_coords);

void quad4_shape_derivatives(double* dN, const double* natural_coords);

void quad4_gauss_points(double* gp_w);

ElementType create_quad4_element();

#endif