#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <direct.h>   /* _getcwd */
#include "mesh.h"
#include "elements.h"
#include "quad4.h"
#include "elasticity_2D_plane_stress.h"
#include "assembly.h"
#include "unv_reader.h"
#include "bc.h"
#include "solver.h"
#include "out_writter.h"
#include "io_plots.h"

/* Tipos de solver disponiveis para a secao de configuracao, dentro de main() */
typedef enum { SOLVER_DENSE, SOLVER_SPARSE } SolverType;

/* Especificacao de uma condicao de contorno de deslocamento: aplica "value"
   no grau de liberdade "dof" (0 = x, 1 = y) de todos os nos do grupo
   "group_name" (grupo definido no pre-processador, ex: SALOME). */
typedef struct { const char* group_name; int dof; double value; } DisplacementBCSpec;

/* Especificacao de uma condicao de contorno de forca: distribui
   "total_force" igualmente entre os nos do grupo "group_name", no grau de
   liberdade "dof" (0 = x, 1 = y). */
typedef struct { const char* group_name; int dof; double total_force; } ForceBCSpec;

int main(void)
{
    /* ================================================================
       SECAO DE CONFIGURACAO
       ================================================================ */

    /* Malha de entrada: arquivo UNV (I-DEAS Universal File) exportado do
       SALOME */
    const char* unv_filename = "C:/Users/lucas/Downloads/Mesh_1.unv";

    /* Solver usado para resolver o sistema global: SOLVER_DENSE (LAPACK
       dsysv) ou SOLVER_SPARSE (MKL PARDISO, direto e paralelo). */
    SolverType solver_type = SOLVER_SPARSE;

    /* Elemento */
    ElementType etype = create_quad4_element();

    /* Fisica */
    PhysicsModel physics;
    physics.dof_per_node = 2;
    physics.integrate = elasticity_2d_plane_stress_integrate;

    /* Propriedades do material (estado plano de tensoes) */
    double E = 200e9;
    double nu = 0.3;
    double material_properties[2] = { E, nu };

    /* Condicoes de contorno de deslocamento, por grupo de nos.
       Os nomes de grupo devem existir na malha UNV (definidos no SALOME). */
    DisplacementBCSpec disp_bc_specs[] = {
        { "BC", 0, 0.0E0 },
        { "BC", 1, 0.0E0 },
    };
    int n_disp_bc_specs = (int)(sizeof(disp_bc_specs) / sizeof(disp_bc_specs[0]));

    /* Condicoes de contorno de forca, por grupo de nos. */
    ForceBCSpec force_bc_specs[] = {
        { "FORCE", 1, 1.0E4 },
    };
    int n_force_bc_specs = (int)(sizeof(force_bc_specs) / sizeof(force_bc_specs[0]));

    /* ================================================================ */

    printf("=== FEM_LIB: malha 2D quad4 (UNV) ===\n\n");

    /* 1. Malha */
    Mesh mesh;
    mesh_init(&mesh);
    if (read_unv_mesh(unv_filename, &mesh) != 0)
    {
        printf("Erro ao ler a malha UNV (%s).\n", unv_filename);
        return 1;
    }
    mesh_print_info(&mesh);

    /* Exporta a malha (sem resultados) para visualizacao imediata */
    write_vtk("malha.vtk", &mesh, NULL);
    printf("Malha exportada para malha.vtk\n");

    int n = mesh.total_dofs;

    printf("\n--- Aplicacao das condicoes de contorno (por grupo) ---\n");

    /* Condicoes de contorno de deslocamento */
    int bc_count = 0;
    for (int s = 0; s < n_disp_bc_specs; s++)
    {
        MeshGroup* g = mesh_find_group(&mesh, disp_bc_specs[s].group_name);
        if (!g)
        {
            printf("Aviso: grupo de deslocamento \"%s\" nao encontrado na malha.\n",
                disp_bc_specs[s].group_name);
            continue;
        }
        //Conta o numero de nos totais com condicoes de contorno de deslocamento
        bc_count += g->n_nodes;
    }

    int* bc_dofs = malloc(bc_count * sizeof(int));
    double* bc_vals = malloc(bc_count * sizeof(double));
    int bc_idx = 0;

    for (int s = 0; s < n_disp_bc_specs; s++)
    {
        MeshGroup* g = mesh_find_group(&mesh, disp_bc_specs[s].group_name);
        if (!g) continue;

        for (int j = 0; j < g->n_nodes; j++)
        {
            //Determina dof global do no "j" do grupo "g" 
            bc_dofs[bc_idx] = g->node_ids[j] * physics.dof_per_node + disp_bc_specs[s].dof;
			//Define o valor do deslocamento prescrito para este dof global
            bc_vals[bc_idx] = disp_bc_specs[s].value;
            bc_idx++;
        }
        printf("Deslocamento: grupo \"%s\", dof %d = %e (%d nos).\n",
            disp_bc_specs[s].group_name, disp_bc_specs[s].dof,
            disp_bc_specs[s].value, g->n_nodes);
    }
	//Cria a variavel disp_bc do tipo DisplacementBC, que armazena o numero de dofs, os indices dos dofs e os valores dos deslocamentos prescritos para ser usado nas rotinas de aplicavao de BC
    DisplacementBC disp_bc = { .n_dofs = bc_idx, .dof = bc_dofs, .value = bc_vals };

    /* Condicoes de contorno de forca */
    int fc_count = 0;
    for (int s = 0; s < n_force_bc_specs; s++)
    {
        MeshGroup* g = mesh_find_group(&mesh, force_bc_specs[s].group_name);
        if (!g)
        {
            printf("Aviso: grupo de forca \"%s\" nao encontrado na malha.\n",
                force_bc_specs[s].group_name);
            continue;
        }
        fc_count += g->n_nodes;
    }

    int* fc_dofs = malloc(fc_count * sizeof(int));
    double* fc_vals = malloc(fc_count * sizeof(double));
    int fc_idx = 0;

    for (int s = 0; s < n_force_bc_specs; s++)
    {
        MeshGroup* g = mesh_find_group(&mesh, force_bc_specs[s].group_name);
        if (!g) continue;

        double per_node = force_bc_specs[s].total_force / g->n_nodes;
        for (int j = 0; j < g->n_nodes; j++)
        {
            fc_dofs[fc_idx] = g->node_ids[j] * physics.dof_per_node + force_bc_specs[s].dof;
            fc_vals[fc_idx] = per_node;
            fc_idx++;
        }
        printf("Forca: grupo \"%s\", dof %d, total %.2e N (%d nos).\n",
            force_bc_specs[s].group_name, force_bc_specs[s].dof,
            force_bc_specs[s].total_force, g->n_nodes);
    }

    ForceBC force_bc = { .n_dofs = fc_idx, .dof = fc_dofs, .value = fc_vals };

    /* Montagem e solucao do sistema global, de acordo com o solver escolhido */
    double* K_global = NULL;
    double* R_global = NULL;
    int* rowIndex = NULL, * columns = NULL;
    double* values = NULL;
    int nnz = 0;

    double* u = calloc(n, sizeof(double));

    if (solver_type == SOLVER_DENSE)
    {
        printf("\n--- Solucao (LAPACK dsysv, denso) ---\n");

        K_global = calloc(n * n, sizeof(double));
        R_global = calloc(n, sizeof(double));

        assemble_global_stiffness(K_global, R_global, &mesh, &etype, &physics, material_properties, NULL);
        apply_displacement_bc_dense(K_global, R_global, n, &disp_bc);
        apply_force_bc(R_global, &force_bc);

        dsolve_sym(n, 1, K_global, R_global);
        memcpy(u, R_global, n * sizeof(double));
    }
    else /* SOLVER_SPARSE */
    {
        printf("\n--- Solucao (MKL PARDISO, esparso) ---\n");

        R_global = calloc(n, sizeof(double));
        apply_force_bc(R_global, &force_bc);

        assemble_global_stiffness_sparse(&rowIndex, &columns, &values, &nnz,
            &mesh, &etype, &physics, material_properties);
        apply_displacement_bc_csr(rowIndex, columns, values, R_global, n, &disp_bc);

        int error = solve_csr_pardiso(n, rowIndex, columns, values, R_global, u, 1);

        if (error == 0)
            printf("PARDISO: solucao obtida com sucesso.\n");
        else
            printf("PARDISO: falha na solucao (codigo de erro = %d).\n", error);
    }

    printf("u[0] (no 0, x) = %e\n", u[0]);
    printf("u[1] (no 0, y) = %e\n", u[1]);
    printf("u[%d] (top-right, y) = %e\n", n - 1, u[n - 1]);

    /* 7. Exportacao VTK (resultados com deslocamento) */
    write_vtk("resultado.vtk", &mesh, u);
    printf("\nResultado exportado para resultado.vtk\n");

    /* 8. Fator de escala do Warp: amplifica a deformada para ~10% do dominio.
       O tamanho do dominio e estimado pela caixa delimitadora (bounding box)
       da malha, ja que ela agora vem de um arquivo externo (UNV). */
    double min_x = mesh.coords[0], max_x = mesh.coords[0];
    double min_y = mesh.coords[1], max_y = mesh.coords[1];
    for (int i = 1; i < mesh.n_nodes; i++)
    {
        double x = mesh.coords[i * 3], y = mesh.coords[i * 3 + 1];
        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
    }
    double domain_size = fmin(max_x - min_x, max_y - min_y);
    double warp_scale = compute_warp_scale(&mesh, u, domain_size);

    /* 9. Gera o script do ParaView com caminhos absolutos */
    char cwd[512];
    if (!_getcwd(cwd, sizeof(cwd)))
        cwd[0] = '\0';

    char mesh_path[600], result_path[600], script_path[600], script_arg[640];
    snprintf(mesh_path, sizeof(mesh_path), "%s\\malha.vtk", cwd);
    snprintf(result_path, sizeof(result_path), "%s\\resultado.vtk", cwd);
    snprintf(script_path, sizeof(script_path), "%s\\view.py", cwd);

    write_paraview_script(script_path, mesh_path, result_path, warp_scale);
    snprintf(script_arg, sizeof(script_arg), "--script=\"%s\"", script_path);

    /* 10. Abre o ParaView ja com a visualizacao montada */
    if (open_in_paraview(script_arg))
        printf("ParaView iniciado (deformada + malha + arestas dos elementos).\n");
    else
        printf("ParaView nao encontrado.\n"
               "  Instale em C:\\Program Files\\ParaView\\\n"
               "  Ou adicione paraview ao PATH e recompila.\n");

    /* 11. Cleanup */
    free(bc_dofs); free(bc_vals);
    free(fc_dofs); free(fc_vals);
    free(u);
    free(K_global); free(R_global);
    free(rowIndex); free(columns); free(values);
    mesh_free(&mesh);

    printf("\n=== FIM ===\n");
    return 0;
}
