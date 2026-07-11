#ifndef OUT_WRITTER_H
#define OUT_WRITTER_H

#include "mesh.h"

/* Escreve a malha (e, opcionalmente, o campo de deslocamento u) em um
   arquivo VTK legacy ASCII (DATASET UNSTRUCTURED_GRID), pronto para abrir
   no ParaView. Se u for NULL, exporta somente a geometria da malha. */
void write_vtk(const char* filename, Mesh* mesh, double* u);

#endif
