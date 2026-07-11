from paraview.simple import *

view = GetActiveViewOrCreate('RenderView')

# Malha de referencia (contorno cinza, nao deformada)
malha = LegacyVTKReader(FileNames=[r'C:\Users\lucas\source\repos\FEM_LIB\FEM_LIB\malha.vtk'])
m = Show(malha, view)
m.Representation = 'Wireframe'
m.AmbientColor = [0.6, 0.6, 0.6]
m.DiffuseColor = [0.6, 0.6, 0.6]

# Resultados: deformada com arestas dos elementos
res = LegacyVTKReader(FileNames=[r'C:\Users\lucas\source\repos\FEM_LIB\FEM_LIB\resultado.vtk'])
warp = WarpByVector(Input=res)
warp.Vectors = ['POINTS', 'displacement']
warp.ScaleFactor = 4975.96
w = Show(warp, view)
w.Representation = 'Surface With Edges'
ColorBy(w, ('POINTS', 'displacement', 'Magnitude'))
w.RescaleTransferFunctionToDataRange(True, False)
w.SetScalarBarVisibility(view, True)

ResetCamera(view)
Render()
