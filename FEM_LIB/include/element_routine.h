#ifndef ELEMENTOROUTINE_H
#define ELEMENTOROUTINE_H

#include "elements.h"
#include "mesh.h"
#include <stdlib.h>
#include <math.h>

void integrate_stiffness_matrix(double* Ke, double* Re, ElementType* element, PhysicsModel* physics, Mesh* mesh, int e, const double* material_properties, const double* u);

#endif