#include "unv_reader.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define UNV_MAX_LINE 512
#define UNV_QUAD4_DESCRIPTOR 44
#define UNV_NODE_ENTITY_TYPE 7

/* Remove espacos/quebras de linha das extremidades da string, em memoria. */
static char* trim(char* s)
{
    while (isspace((unsigned char)*s)) s++;
    if (*s == '\0') return s;

    char* end = s + strlen(s) - 1;
    while (end > s && isspace((unsigned char)*end)) *end-- = '\0';

    return s;
}

/* Converte o expoente Fortran ('1.0D+00') para notacao C ('1.0E+00'), pois
   varios exportadores UNV usam 'D' e sscanf nao reconhece esse expoente. */
static void fortran_to_c_exponent(char* line)
{
    for (int i = 0; line[i]; i++)
        if (line[i] == 'D' || line[i] == 'd')
            line[i] = 'E';
}

/* Vetores dinamicos auxiliares (capacidade dobra conforme necessario),
   usados pois o numero de nos/elementos/entidades de grupo so e conhecido
   ao final da leitura de cada bloco. */
typedef struct { double* data; int size, cap; } DVec;
typedef struct { int* data; int size, cap; } IVec;

static void dvec_push3(DVec* v, double a, double b, double c)
{
    if (v->size + 3 > v->cap)
    {
        v->cap = v->cap ? v->cap * 2 : 256;
        v->data = realloc(v->data, v->cap * sizeof(double));
    }
    v->data[v->size++] = a;
    v->data[v->size++] = b;
    v->data[v->size++] = c;
}

static void ivec_push(IVec* v, int x)
{
    if (v->size + 1 > v->cap)
    {
        v->cap = v->cap ? v->cap * 2 : 256;
        v->data = realloc(v->data, v->cap * sizeof(int));
    }
    v->data[v->size++] = x;
}

/* Le todos os inteiros de uma linha para "out", retornando a quantidade lida. */
static int parse_ints(const char* line, int* out, int max_out)
{
    int count = 0;
    const char* p = line;

    while (*p && count < max_out)
    {
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0' || *p == '\n' || *p == '\r') break;

        int val, consumed;
        if (sscanf_s(p, "%d%n", &val, &consumed) != 1) break;

        out[count++] = val;
        p += consumed;
    }

    return count;
}

/* Bloco de nos (dataset 2411). Cada no ocupa 2 linhas: uma linha de
   cabecalho (label, sistema de coordenadas, cor - ignorados aqui) e uma
   linha com x, y, z. Termina na linha "-1". */
static void read_nodes_block(FILE* fp, DVec* coords)
{
    char line[UNV_MAX_LINE];

    while (fgets(line, sizeof(line), fp))
    {
        char* t = trim(line);
        if (strcmp(t, "-1") == 0) return;

        if (!fgets(line, sizeof(line), fp)) return;

        fortran_to_c_exponent(line);
        double x = 0.0, y = 0.0, z = 0.0;
        sscanf_s(line, "%lf %lf %lf", &x, &y, &z);

        dvec_push3(coords, x, y, z);
    }
}

/* Bloco de elementos (dataset 2412). Cada elemento ocupa 2 linhas:
   cabecalho (label, fe_descriptor, phys/mat table, cor, num_nos) e a
   lista de nos. Mantem apenas elementos QUAD4 (fe_descriptor == 44);
   os demais (vigas, linhas de grupo, etc.) sao descartados. */
static void read_elements_block(FILE* fp, IVec* conn)
{
    char line[UNV_MAX_LINE];

    while (fgets(line, sizeof(line), fp))
    {
        char* t = trim(line);
        if (strcmp(t, "-1") == 0) return;

        int label, fe_descriptor, phys_tab, mat_tab, color, num_nodes;
        int n_read = sscanf_s(line, "%d %d %d %d %d %d",
            &label, &fe_descriptor, &phys_tab, &mat_tab, &color, &num_nodes);

        if (n_read != 6) continue;
        if (!fgets(line, sizeof(line), fp)) return;

        if (fe_descriptor == UNV_QUAD4_DESCRIPTOR && num_nodes == 4)
        {
            int n0, n1, n2, n3;
            sscanf_s(line, "%d %d %d %d", &n0, &n1, &n2, &n3);

            /* UNV usa indices de no 1-based; a malha usa 0-based */
            ivec_push(conn, n0 - 1);
            ivec_push(conn, n1 - 1);
            ivec_push(conn, n2 - 1);
            ivec_push(conn, n3 - 1);
        }
    }
}

/* Bloco de grupos permanentes (dataset 2467/2477). Cada grupo tem um
   cabecalho de 8 inteiros (o ultimo e o numero de entidades), uma linha
   de nome, e quadruplas (tipo, tag, node_leaf_id, componente) de
   entidades - 2 entidades por linha (8 inteiros). Mantem somente
   entidades do tipo "no" (codigo 7); entidades de outro tipo (elementos,
   etc.) sao contadas para consumir a quantidade correta de registros,
   mas descartadas. */
static void read_groups_block(FILE* fp, MeshGroup** groups_out, int* n_groups_out)
{
    char line[UNV_MAX_LINE];
    int cap = 8, n_groups = 0;
    MeshGroup* groups = malloc(cap * sizeof(MeshGroup));

    while (fgets(line, sizeof(line), fp))
    {
        char* t = trim(line);
        if (strcmp(t, "-1") == 0) break;

        int header[8];
        if (parse_ints(line, header, 8) < 8) break;
        int num_entities = header[7];

        char name_line[UNV_MAX_LINE];
        if (!fgets(name_line, sizeof(name_line), fp)) break;
        char* name = trim(name_line);

        IVec node_ids = { 0 };
        int collected = 0;
        while (collected < num_entities && fgets(line, sizeof(line), fp))
        {
            int ints[64];
            int n_ints = parse_ints(line, ints, 64);

            for (int i = 0; i + 3 < n_ints && collected < num_entities; i += 4, collected++)
            {
                int type = ints[i];
                int tag = ints[i + 1];
                if (type == UNV_NODE_ENTITY_TYPE)
                    ivec_push(&node_ids, tag - 1);
            }
        }

        if (n_groups == cap)
        {
            cap *= 2;
            groups = realloc(groups, cap * sizeof(MeshGroup));
        }

        strncpy_s(groups[n_groups].name, sizeof(groups[n_groups].name), name, _TRUNCATE);
        groups[n_groups].n_nodes = node_ids.size;
        groups[n_groups].node_ids = node_ids.data;
        n_groups++;
    }

    *groups_out = groups;
    *n_groups_out = n_groups;
}

int read_unv_mesh(const char* filename, Mesh* mesh)
{
    FILE* fp = fopen(filename, "r");
    if (!fp)
    {
        printf("Cannot open file\n");
        return -1;
    }

    DVec coords = { 0 };
    IVec conn = { 0 };
    MeshGroup* groups = NULL;
    int n_groups = 0;

    char line[UNV_MAX_LINE];
    while (fgets(line, sizeof(line), fp))
    {
        char* t = trim(line);
        if (strcmp(t, "-1") != 0)
            continue;

        if (!fgets(line, sizeof(line), fp)) break;
        int dataset_id = atoi(trim(line));

        if (dataset_id == 2411)
            read_nodes_block(fp, &coords);
        else if (dataset_id == 2412)
            read_elements_block(fp, &conn);
        else if (dataset_id == 2467 || dataset_id == 2477)
            read_groups_block(fp, &groups, &n_groups);
        else
        {
            /* bloco desconhecido/nao suportado: apenas consome ate o "-1" de fechamento */
            while (fgets(line, sizeof(line), fp))
                if (strcmp(trim(line), "-1") == 0) break;
        }
    }

    fclose(fp);

    if (coords.size == 0 || conn.size == 0)
    {
        free(coords.data);
        free(conn.data);
        for (int i = 0; i < n_groups; i++) free(groups[i].node_ids);
        free(groups);
        printf("UNV: nenhum no ou elemento QUAD4 encontrado\n");
        return -1;
    }

    mesh->dim = 2;
    mesh->n_nodes = coords.size / 3;
    mesh->n_elements = conn.size / 4;
    mesh->nodes_per_element = 4;
    mesh->dof_per_node = 2;
    mesh->total_dofs = mesh->n_nodes * mesh->dof_per_node;
    mesh->coords = coords.data;
    mesh->connectivity = conn.data;
    mesh->groups = groups;
    mesh->n_groups = n_groups;

    return 0;
}
