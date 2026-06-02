#ifndef ELEMENTOROUTINE_H
#define ELEMENTOROUTINE_H

#include "elements.h"
#include <stdlib.h>
#include <math.h>

void integrate_stiffness_matrix(double* Ke, ElementType* element, const double* node_coords, const double* material_properties);

#endif