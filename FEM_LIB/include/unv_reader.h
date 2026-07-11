#ifndef UNV_READER_H
#define UNV_READER_H

#include "mesh.h"

/* Le um arquivo UNV (I-DEAS Universal File), tipicamente exportado pelo
   SALOME, e preenche a malha: nos (dataset 2411), elementos quad4
   (dataset 2412, fe descriptor 44) e grupos de nos (dataset 2467/2477),
   usados para aplicar condicoes de contorno por nome.
   Elementos que nao sejam quad4 sao ignorados.
   Retorna 0 em sucesso, -1 em erro. */
int read_unv_mesh(const char* filename, Mesh* mesh);

#endif
