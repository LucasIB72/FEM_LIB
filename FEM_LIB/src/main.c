#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <io.h>
#include "mesh.h"
#include "elements.h"
#include "quad4.h"
#include "elasticity_2D_plane_stress.h"
#include "assembly.h"
#include "simp_input.h"
#include "bc.h"
#include "solver.h"

static void write_vtk(const char* filename, Mesh* mesh, double* u)
{
    FILE* fp = fopen(filename, "w");
    if (!fp) return;

    fprintf(fp, "# vtk DataFile Version 3.0\nFEM_LIB\nASCII\nDATASET UNSTRUCTURED_GRID\n");
    fprintf(fp, "POINTS %d double\n", mesh->n_nodes);
    for (int i = 0; i < mesh->n_nodes; i++)
        fprintf(fp, "%e %e %e\n", mesh->coords[i * 3], mesh->coords[i * 3 + 1], mesh->coords[i * 3 + 2]);

    fprintf(fp, "CELLS %d %d\n", mesh->n_elements, mesh->n_elements * (mesh->nodes_per_element + 1));
    for (int e = 0; e < mesh->n_elements; e++)
    {
        fprintf(fp, "%d", mesh->nodes_per_element);
        for (int j = 0; j < mesh->nodes_per_element; j++)
            fprintf(fp, " %d", mesh->connectivity[e * mesh->nodes_per_element + j]);
        fprintf(fp, "\n");
    }

    fprintf(fp, "CELL_TYPES %d\n", mesh->n_elements);
    for (int e = 0; e < mesh->n_elements; e++)
        fprintf(fp, "%d\n", mesh->nodes_per_element == 4 ? 9 : 5);

    if (u)
    {
        fprintf(fp, "POINT_DATA %d\n", mesh->n_nodes);
        fprintf(fp, "VECTORS displacement double\n");
        for (int i = 0; i < mesh->n_nodes; i++)
            fprintf(fp, "%e %e 0.0\n", u[i * 2], u[i * 2 + 1]);
    }

    fclose(fp);
}

int main(void)
{
    printf("=== FEM_LIB: malha 2D quad4 ===\n\n");

    double L_x = 1.0, L_y = 1.0, elem_size = 0.1;

    /* 1. Malha */
    Mesh mesh;
    mesh_init(&mesh);
    if (generate_rect_mesh_to_mesh(L_x, L_y, elem_size, &mesh) != 0)
    {
        printf("Erro ao gerar a malha.\n");
        return 1;
    }
    mesh_print_info(&mesh);

    /* 2. Elemento */
    ElementType etype = create_quad4_element();

    /* 3. Fisica (plane stress) */
    PhysicsModel physics;
    physics.dof_per_node = 2;
    physics.integrate = elasticity_2d_plane_stress_integrate;

    /* 4. Propriedades do material */
    double material_properties[2] = { 200e9, 0.3 };

    /* 5. Alocacao dos sistemas globais */
    int n = mesh.total_dofs;
    double* K_global = calloc(n * n, sizeof(double));
    double* R_global = calloc(n, sizeof(double));

    if (!K_global || !R_global)
    {
        printf("Falha ao alocar matrizes globais.\n");
        free(K_global); free(R_global);
        mesh_free(&mesh);
        return 1;
    }

    /* 6. Montagem */
    assemble_global_stiffness(K_global, R_global, &mesh, &etype, &physics, material_properties, NULL);

    printf("\n--- Aplicacao das condicoes de contorno ---\n");

    /* 7. Condicoes de contorno de deslocamento (aresta esquerda: x = 0) */
    int* bc_dofs = malloc(mesh.n_nodes * 2 * sizeof(int));
    double* bc_vals = malloc(mesh.n_nodes * 2 * sizeof(double));
    int bc_count = 0;

    for (int i = 0; i < mesh.n_nodes; i++)
    {
        if (fabs(mesh.coords[i * 3]) < 1e-12)
        {
            bc_dofs[bc_count] = i * 2;
            bc_vals[bc_count] = 0.0;
            bc_count++;
            bc_dofs[bc_count] = i * 2 + 1;
            bc_vals[bc_count] = 0.0;
            bc_count++;
        }
    }

    DisplacementBC disp_bc = { .n_dofs = bc_count, .dof = bc_dofs, .value = bc_vals };
    apply_displacement_bc_dense(K_global, R_global, n, &disp_bc);
    printf("Deslocamento: %d graus de liberdade prescritos (x = 0).\n", bc_count);

    /* 8. Condicoes de contorno de forca (aresta superior: y = L_y) */
    int n_bc_f = 0;
    for (int i = 0; i < mesh.n_nodes; i++)
        if (fabs(mesh.coords[i * 3 + 1] - L_y) < 1e-12)
            n_bc_f++;

    int* fc_dofs = malloc(n_bc_f * sizeof(int));
    double* fc_vals = malloc(n_bc_f * sizeof(double));
    double F_total = 1e6;
    int fc_count = 0;

    for (int i = 0; i < mesh.n_nodes; i++)
    {
        if (fabs(mesh.coords[i * 3 + 1] - L_y) < 1e-12)
        {
            fc_dofs[fc_count] = i * 2 + 1;
            fc_vals[fc_count] = F_total / n_bc_f;
            fc_count++;
        }
    }

    ForceBC force_bc = { .n_dofs = n_bc_f, .dof = fc_dofs, .value = fc_vals };
    apply_force_bc(R_global, &force_bc);
    printf("Forca: %d nos com carga vertical total de %.2e N.\n", n_bc_f, F_total);

    /* 9. Solucao do sistema denso */
    printf("\n--- Solucao (dense LU) ---\n");
    double* u = calloc(n, sizeof(double));
    double* K_copy = malloc(n * n * sizeof(double));
    memcpy(K_copy, K_global, n * n * sizeof(double));

    dsolve_sym(n, 1, K_copy, R_global);
    memcpy(u, R_global, n * sizeof(double));
    free(K_copy);

    printf("u[0] (no 0, x) = %e\n", u[0]);
    printf("u[1] (no 0, y) = %e\n", u[1]);
    printf("u[%d] (top-right, y) = %e\n", n - 1, u[n - 1]);

    /* 10. Solucao do sistema esparso com CG */
    printf("\n--- Solucao (sparse CG) ---\n");

    int* rowIndex = NULL, *columns = NULL;
    double* values = NULL;
    int nnz = 0;

    double* F_sparse = calloc(n, sizeof(double));
    double* K_sparse = calloc(n * n, sizeof(double));

    assemble_global_stiffness(K_sparse, F_sparse, &mesh, &etype, &physics, material_properties, NULL);
    apply_displacement_bc_dense(K_sparse, F_sparse, n, &disp_bc);
    apply_force_bc(F_sparse, &force_bc);

    assemble_global_stiffness_sparse(&rowIndex, &columns, &values, &nnz,
        &mesh, &etype, &physics, material_properties);
    apply_displacement_bc_csr(rowIndex, columns, values, F_sparse, n, &disp_bc);

    double* u_cg = calloc(n, sizeof(double));
    int iter = solve_csr_cg(n, rowIndex, columns, values, F_sparse, u_cg, 2000, 1e-12);

    if (iter >= 0)
        printf("CG convergiu em %d iteracoes.\n", iter);
    else
        printf("CG nao convergiu.\n");

    printf("u_cg[0] (no 0, x) = %e\n", u_cg[0]);
    printf("u_cg[1] (no 0, y) = %e\n", u_cg[1]);
    printf("u_cg[%d] (top-right, y) = %e\n", n - 1, u_cg[n - 1]);

    printf("\nDiferenca max LU vs CG: %e\n",
        fabs(u[n - 1] - u_cg[n - 1]));

    /* 11. Exportacao VTK */
    write_vtk("resultado.vtk", &mesh, u);
    printf("\nResultado exportado para resultado.vtk\n");

    /* 12. Abre ParaView automaticamente com visualizacao */
    {
        printf("Abrindo ParaView...\n");
        char cmd_pv[1024];
        int launched = 0;

        const char* pv_paths[] = {
            "C:\\Program Files\\ParaView\\bin\\paraview.exe",
            "C:\\Program Files\\ParaView 6.1.1\\bin\\paraview.exe",
            "C:\\Program Files\\ParaView 5.12.0\\bin\\paraview.exe",
            "C:\\Program Files\\ParaView 5.11.0\\bin\\paraview.exe",
            "C:\\Program Files\\ParaView 5.10.0\\bin\\paraview.exe",
        };
        for (int i = 0; i < 5 && !launched; i++)
        {
            if (_access(pv_paths[i], 0) == 0)
            {
                sprintf(cmd_pv, "\"%s\" --script=auxiliary\\plot_results.py resultado.vtk",
                    pv_paths[i]);
                system(cmd_pv);
                launched = 1;
            }
        }

        if (launched)
            printf("ParaView iniciado.\n");
        else
            printf("ParaView nao encontrado.\n  Instale em C:\\Program Files\\ParaView\\\n"
                   "  Ou abra manualmente:\n"
                   "    paraview resultado.vtk\n");
    }

    /* 13. Cleanup */
    free(bc_dofs); free(bc_vals);
    free(fc_dofs); free(fc_vals);
    free(u); free(u_cg);
    free(K_global); free(R_global);
    free(K_sparse); free(F_sparse);
    free(rowIndex); free(columns); free(values);
    mesh_free(&mesh);

    printf("\n=== FIM ===\n");
    return 0;
}
