#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inp_reader.h"
#include "mesh.h"

#define MAX_LINE 1024

int read_inp_mesh(const char* filename,
    double** coords,
    int** connectivity,
    int* n_nodes,
    int* n_elements,
    int* nodes_per_element)
{
    FILE* fp = fopen(filename, "r");

    if (!fp)
    {
        printf("Cannot open file\n");
        return -1;
    }

    char line[MAX_LINE];

    int reading_nodes = 0;
    int reading_elements = 0;

    int node_count = 0;
    int elem_count = 0;
    int elem_node_count = 0;

	int node_index = 0;
	int elem_index = 0;

	// Contagem do número de nós e elementos para alocar memória posteriormente
    while (fgets(line, MAX_LINE, fp))
    {

        if ((strncmp(line, "*Node", 5) == 0) || (strncmp(line, "*NODE", 5) == 0))
        {
            reading_nodes = 1;
            reading_elements = 0;
            continue;
        }

        if ((strncmp(line, "*Element", 8) == 0) || (strncmp(line, "*ELEMENT", 8) == 0))
        {
            reading_nodes = 0;
            reading_elements = 1;
            continue;
        }

        if (line[0] == '*')
        {
            reading_nodes = 0;
            reading_elements = 0;
            continue;
        }

        if (reading_nodes)
            node_count++;

        if (reading_elements)
        {
            if (elem_count == 0)
            {
                int values[32];

                for (int i = 0; line[i]; i++)
                    if (line[i] == ',')
                        line[i] = ' ';

                char* ptr = line;

                // Lê os números inteiros "%d" da linha e armazena em values[]
                while (sscanf_s(ptr, "%d", &values[elem_node_count]) == 1)
                {
                    // Move o ponteiro para o final do número lido
                    while (*ptr && *ptr != ' ')
                        ptr++;
                    // Faz o ponteiro pular os espaços em branco
                    while (*ptr == ' ')
                        ptr++;
					// Verifica se chegamos ao final da linha
					if (*ptr == '\0')
                        break;
                    elem_node_count++;
                }
                elem_node_count = elem_node_count - 1;
            }
            elem_count++;
        }
    }

    rewind(fp);

    reading_nodes = 0;
    reading_elements = 0;

    *n_nodes = node_count;
    *n_elements = elem_count;
    *nodes_per_element = elem_node_count;

	// Salvando matriz de coordenadas e conectividade
    *coords = malloc((*n_nodes) * 3 * sizeof(double));
    *connectivity = malloc((*n_elements) * (*nodes_per_element) * sizeof(int));

    while (fgets(line, MAX_LINE, fp))
    {
        if ((strncmp(line, "*Node", 5) == 0) || (strncmp(line, "*NODE", 5) == 0))
        {
            reading_nodes = 1;
            reading_elements = 0;
            continue;
        }

        if ((strncmp(line, "*Element", 8) == 0) || (strncmp(line, "*ELEMENT", 8) == 0))
        {
            reading_nodes = 0;
            reading_elements = 1;
            continue;
        }

        if (line[0] == '*')
        {
            reading_nodes = 0;
            reading_elements = 0;
            continue;
        }

        // Armazena as coordenadas dos nós dos elementos no vetor coords

        if (reading_nodes)
        {
            int id;
            double x, y, z;

            for (int i = 0; line[i]; i++)
                if (line[i] == ',')
                    line[i] = ' ';

            sscanf_s(line, "%d %lf %lf %lf", &id, &x, &y, &z);

            (*coords)[node_index * 3 + 0] = x;
            (*coords)[node_index * 3 + 1] = y;
            (*coords)[node_index * 3 + 2] = z;

            node_index++;
        }
        
		// Armazena a conectividade dos elementos no vetor connectivity

        if (reading_elements)
        {
            int values[32];
            int count = 0;

            for (int i = 0; line[i]; i++)
                if (line[i] == ',')
                    line[i] = ' ';

            char* ptr = line;

            while (sscanf_s(ptr, "%d", &values[count]) == 1)
            {
                while (*ptr && *ptr != ' ')
                    ptr++;

                while (*ptr == ' ')
                    ptr++;

                if (*ptr == '\0')
                    break;

                count++;
            }

            /* first value = element id */
            for (int i = 0; i < elem_node_count; i++)
            {
                (*connectivity)[elem_index * elem_node_count + i] = values[i + 1] - 1;
            }

            elem_index++;
        }
    }

    fclose(fp);

    return 0;
}