#include "io_plots.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double compute_warp_scale(const Mesh* mesh, const double* u, double domain_size)
{
    double max_mag = 0.0;
    for (int i = 0; i < mesh->n_nodes; i++)
    {
        double ux = u[i * 2], uy = u[i * 2 + 1];
        double mag = sqrt(ux * ux + uy * uy);
        if (mag > max_mag) max_mag = mag;
    }

    if (max_mag <= 0.0) return 1.0;

    return 0.1 * domain_size / max_mag;
}

void write_paraview_script(const char* script_path,
                            const char* mesh_vtk,
                            const char* result_vtk,
                            double warp_scale)
{
    FILE* fp = fopen(script_path, "w");
    if (!fp) return;

    fprintf(fp,
        "from paraview.simple import *\n"
        "\n"
        "view = GetActiveViewOrCreate('RenderView')\n"
        "\n"
        "# Malha de referencia (contorno cinza, nao deformada)\n"
        "malha = LegacyVTKReader(FileNames=[r'%s'])\n"
        "m = Show(malha, view)\n"
        "m.Representation = 'Wireframe'\n"
        "m.AmbientColor = [0.6, 0.6, 0.6]\n"
        "m.DiffuseColor = [0.6, 0.6, 0.6]\n"
        "\n"
        "# Resultados: deformada com arestas dos elementos\n"
        "res = LegacyVTKReader(FileNames=[r'%s'])\n"
        "warp = WarpByVector(Input=res)\n"
        "warp.Vectors = ['POINTS', 'displacement']\n"
        "warp.ScaleFactor = %g\n"
        "w = Show(warp, view)\n"
        "w.Representation = 'Surface With Edges'\n"
        "ColorBy(w, ('POINTS', 'displacement', 'Magnitude'))\n"
        "w.RescaleTransferFunctionToDataRange(True, False)\n"
        "w.SetScalarBarVisibility(view, True)\n"
        "\n"
        "ResetCamera(view)\n"
        "Render()\n",
        mesh_vtk, result_vtk, warp_scale);

    fclose(fp);
}

int open_in_paraview(const char* args)
{
    char cmd[1024];

    /* 1. paraview disponivel no PATH */
    snprintf(cmd, sizeof(cmd),
        "where paraview >nul 2>nul && start \"\" /B paraview %s", args);
    if (system(cmd) == 0)
        return 1;

    /* 2. procura qualquer instalacao versionada em Program Files.
          O 'exit /b 0' so e alcancado quando um paraview.exe e encontrado;
          caso contrario o 'exit /b 1' final sinaliza falha. */
    snprintf(cmd, sizeof(cmd),
        "for /d %%D in (\"C:\\Program Files\\ParaView*\") do "
        "@if exist \"%%D\\bin\\paraview.exe\" "
        "( start \"\" /B \"%%D\\bin\\paraview.exe\" %s & exit /b 0 ) "
        "& exit /b 1",
        args);
    if (system(cmd) == 0)
        return 1;

    return 0;
}
