#ifndef ELEMENTOROUTINE_H
#define ELEMENTOROUTINE_H

#include "mesh.h"
#include "elements.h"
#include "materials.h"

typedef struct ElementRoutine
{
    void (*compute)(
        double* Ke,
        double* Re,
        Mesh* mesh,
        ElementType* etype,
        Material* mat,
        int element_id,
        double* u);

} ElementRoutine;

#endif 