# Vista de ParaView para Kuntur.
# Genera paraview/kuntur.pvsm (File → Load State).
# También se puede correr con:
#   pvpython paraview/kuntur_vista.py

from paraview.simple import *
import os

root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
vtu_path = os.path.join(root, "results", "helmholtz_magnetorquers.vtu")
torque_path = os.path.join(root, "results", "torque.vtk")
state_path = os.path.join(root, "paraview", "kuntur.pvsm")

if not os.path.isfile(vtu_path):
    raise SystemExit("Falta " + vtu_path + ". Corre primero julia --project=. main.jl")

# --- Datos ----------------------------------------------------------
campo = XMLUnstructuredGridReader(registrationName="Campo", FileName=[vtu_path])
campo.UpdatePipeline()

torque = LegacyVTKReader(registrationName="Torque", FileNames=[torque_path])
flecha = Tube(registrationName="FlechaTorque", Input=torque)
flecha.Radius = 0.012
flecha.NumberofSides = 12

# Plano XZ: corta las dos bobinas de Helmholtz y las varillas
corte = Slice(registrationName="PlanoXZ", Input=campo)
corte.SliceType = "Plane"
corte.SliceType.Origin = [0.0, 0.0, 0.0]
corte.SliceType.Normal = [0.0, 1.0, 0.0]

# Intensidad a lo largo del eje del par
linea_z = PlotOverLine(registrationName="LineaEjeZ", Input=campo)
linea_z.Point1 = [0.0, 0.0, -1.7]
linea_z.Point2 = [0.0, 0.0, 1.7]
linea_z.SamplingPattern = "Sample Uniformly"
linea_z.Resolution = 800
eje_z = Calculator(registrationName="IntensidadEjeZ", Input=linea_z)
eje_z.Function = "coordsZ"
eje_z.ResultArrayName = "z"

# Intensidad a lo largo de la varilla en X (z = z_mid)
linea_x = PlotOverLine(registrationName="LineaEjeX", Input=campo)
linea_x.Point1 = [-0.80, 0.0, 0.039]
linea_x.Point2 = [0.80, 0.0, 0.039]
linea_x.SamplingPattern = "Sample Uniformly"
linea_x.Resolution = 600
eje_x = Calculator(registrationName="IntensidadEjeX", Input=linea_x)
eje_x.Function = "coordsX"
eje_x.ResultArrayName = "x"

# --- Layout ---------------------------------------------------------
layout = CreateLayout(name="Kuntur")
vista3d = CreateView("RenderView")
grafica_z = CreateView("XYChartView")
grafica_x = CreateView("XYChartView")

layout.AssignView(0, vista3d)
layout.SplitHorizontal(0, 0.58)
layout.AssignView(2, grafica_z)
layout.SplitVertical(2, 0.55)
layout.AssignView(6, grafica_x)
layout.SetSize(1600, 900)

# --- Render ---------------------------------------------------------
corte_disp = Show(corte, vista3d, "GeometryRepresentation")
corte_disp.Representation = "Surface"
ColorBy(corte_disp, ("POINTS", "normB_helmholtz"))
corte_disp.SetScalarBarVisibility(vista3d, True)
corte_disp.RescaleTransferFunctionToDataRange(False, True)

lut = GetColorTransferFunction("normB_helmholtz")
lut.ApplyPreset("Viridis (matplotlib)", True)
bar = GetScalarBar(lut, vista3d)
bar.Title = "B Helmholtz (T)"
bar.ComponentTitle = ""

flecha_disp = Show(flecha, vista3d, "GeometryRepresentation")
flecha_disp.Representation = "Surface"
ColorBy(flecha_disp, None)
flecha_disp.AmbientColor = [1.0, 0.82, 0.15]
flecha_disp.DiffuseColor = [1.0, 0.82, 0.15]

vista3d.CameraPosition = [0.15, -2.4, 0.35]
vista3d.CameraFocalPoint = [0.0, 0.0, 0.0]
vista3d.CameraViewUp = [0.0, 0.0, 1.0]
vista3d.Background = [0.16, 0.17, 0.20]
vista3d.OrientationAxesVisibility = 1

# --- Gráfica |B|(z) -------------------------------------------------
gz = Show(eje_z, grafica_z, "XYChartRepresentation")
gz.XArrayName = "z"
gz.SeriesVisibility = ["normB_helmholtz"]
gz.SeriesColor = ["normB_helmholtz", "0.95", "0.55", "0.15"]
grafica_z.ChartTitle = "Intensidad en el eje del par"
grafica_z.BottomAxisTitle = "z (m)"
grafica_z.LeftAxisTitle = "|B| (T)"
grafica_z.ShowLegend = 1

# --- Gráfica |B|(x) sobre la varilla --------------------------------
gx = Show(eje_x, grafica_x, "XYChartRepresentation")
gx.XArrayName = "x"
gx.SeriesVisibility = ["normB_helmholtz", "normB"]
gx.SeriesColor = [
    "normB_helmholtz", "0.95", "0.55", "0.15",
    "normB", "0.35", "0.65", "0.95",
]
grafica_x.ChartTitle = "Intensidad a lo largo del magnetorquer X"
grafica_x.BottomAxisTitle = "x (m)"
grafica_x.LeftAxisTitle = "|B| (T)"
grafica_x.ShowLegend = 1

SaveState(state_path)
print("Estado guardado en", state_path)
