#pragma once
#include "mesh.h"
#ifndef INP_READER_H
#define INP_READER_H


int read_inp_mesh(const char* filename,
    double** coords,
    int** connectivity,
    int* n_nodes,
    int* n_elements,
    int* nodes_per_element);

#endif