PROJ_NAME := julia-bemholtz


all: build

paraview: build
	cd build && cpack -G TGZ


gmsh: gmsh 


