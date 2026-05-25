
PROJ_NAME := julia-bemholtz

MSH_NAME := torus_coil.geo
MSH_FOLDER := meshes

all: build

paraview: 
	paraview 

gmsh: 
	gmsh ${MSH_FOLDER}/${MSH_NAME}


