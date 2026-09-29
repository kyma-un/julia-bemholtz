
PROJ_NAME := julia-bemholtz

MSH_NAME := torus_coil
MSH_FOLDER := gmsh

all: build

paraview: 
	paraview 

mesh:
	gmsh ${MSH_FOLDER}/${MSH_NAME}.geo

clean_results:
	println("Borrando resultados en /results...")
	rm results/*.vtu

sim:
	julia --project=. main.jl
