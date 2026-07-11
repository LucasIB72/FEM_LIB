#ifndef IO_PLOTS_H
#define IO_PLOTS_H

#include "mesh.h"

/* Calcula um fator de escala de Warp By Vector que amplia o maior
   deslocamento nodal para ~10% do menor lado do dominio (domain_size),
   de forma que a deformada fique visivelmente distinguivel da malha
   original. Retorna 1.0 se nao houver deslocamento (evita escala nula). */
double compute_warp_scale(const Mesh* mesh, const double* u, double domain_size);

/* Gera um script Python do ParaView que carrega a malha e os resultados,
   aplica Warp By Vector (deformada), mostra as arestas dos elementos e
   colore pela magnitude do deslocamento. */
void write_paraview_script(const char* script_path,
                            const char* mesh_vtk,
                            const char* result_vtk,
                            double warp_scale);

/* Abre um ou mais arquivos/argumentos no ParaView (tenta o PATH e depois
   qualquer instalacao encontrada em Program Files). "args" e uma string
   ja pronta (ex: "malha.vtk resultado.vtk" ou "--script=\"view.py\"").
   Retorna 1 se o ParaView foi iniciado, 0 caso contrario. */
int open_in_paraview(const char* args);

#endif
