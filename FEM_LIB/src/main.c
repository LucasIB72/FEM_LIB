#include <stdio.h>
#include "inp_reader.h"
#include "mesh.h"

/*

Relembrando ponteiros :

- & signinifica "endereço de" e é usado para obter o endereço de uma variável.
- * é usado para declarar um ponteiro e para acessar o valor apontado por um ponteiro.
- -> é usado para acessar membros de uma estrutura (apenas estruturas) através de um ponteiro.

Exemplo:
int x = 10;
int a;
int* p = &x; // p é um ponteiro para o endereço de memória de x
printf("%d\n", *p); // Imprime o valor de x (10)
printf("%p\n", p); // Imprime o endereço de memória de x
a = *p; // Atribui o valor apontado por p (10)
*/ 


int main(void)
{
    char filename[] = "C:\\Users\\lucas\\OneDrive\\Desktop\\FEM\\r1.inp";

    double* coords;
    int* connectivity;

    int n_nodes;
    int n_elements;
    int nodes_per_element;

    read_inp_mesh(
        filename,
        &coords,
        &connectivity,
        &n_nodes,
        &n_elements,
        &nodes_per_element
    );

    for (int i = 0; i < n_nodes; i++)
    {
        printf("%d : %f %f %f\n",
            i,
            coords[i * 3 + 0],
            coords[i * 3 + 1],
            coords[i * 3 + 2]);
    }

    return 0;
}
