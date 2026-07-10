from paraview.simple import *
import sys, os

# Abre o VTK
vtk_file = sys.argv[1] if len(sys.argv) > 1 else "resultado.vtk"

mesh = LegacyVTKReader(FileNames=[vtk_file])

# View
view = CreateRenderView()
view.ViewSize = [1200, 800]

# Malha original com edges
disp = Show(mesh, view)
disp.Representation = "Surface With Edges"
ColorBy(disp, ("POINTS", "displacement"))
disp.RescaleTransferFunctionToDataRange(True, False)

# Warp por deslocamento
warp = WarpByVector(Input=mesh)
warp.Vectors = ["POINTS", "displacement"]
warp.ScaleFactor = 1.0

disp_w = Show(warp, view)
disp_w.Representation = "Surface"
ColorBy(disp_w, ("POINTS", "displacement"))
disp_w.RescaleTransferFunctionToDataRange(True, False)

# Ajusta camera
view.ResetCamera()

# Salva screenshot na area de trabalho
desktop = os.path.join(os.path.expanduser("~"), "Desktop")
png_path = os.path.join(desktop, "resultado_FEM.png")
SaveScreenshot(png_path, view, magnification=2)
print(f"Screenshot salvo em {png_path}")

# Abre a interface grafica
interact()
