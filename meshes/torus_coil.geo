// ============================================================
// Bobina de Helmholtz — malla 3D volumétrica para FEM
// ============================================================
SetFactory("OpenCASCADE");

// --- Parámetros ---------------------------------------------
R     = 0.100;    // Radio principal de la bobina [m]
r_w   = 0.008;    // Radio sección transversal del hilo [m]
d     = R / 2;    // Separación centro→bobina (condición Helmholtz: d = R/2)
R_dom = 3.0 * R;  // Radio del dominio de aire

// --- Geometría ----------------------------------------------

// Bobina superior
Torus(1) = {0, 0,  d, R, r_w};

// Bobina inferior
Torus(2) = {0, 0, -d, R, r_w};

// Dominio de aire (esfera)
Sphere(3) = {0, 0, 0, R_dom};

// Fragmentar — crea interfaces coherentes entre bobinas y aire
// Resultado: coil1(1), coil2(2), aire(3)
BooleanFragments{ Volume{3}; Delete; }{ Volume{1}; Volume{2}; Delete; }

// --- Grupos físicos -----------------------------------------
Physical Volume("coil1")     = {1};
Physical Volume("coil2")     = {2};
Physical Volume("air")       = {3};
Physical Surface("boundary") = {1};   // superficie exterior esfera

// --- Control de tamaño de malla -----------------------------
Mesh.Algorithm3D = 4;   // Frontal-Delaunay → tetraedros limpios

// Fino en bobinas, grueso lejos
Field[1] = Distance;
Field[1].VolumesList = {1, 2};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = r_w / 2;        // fino dentro de la bobina
Field[2].SizeMax = R_dom / 5;      // grueso en el exterior
Field[2].DistMin = r_w;
Field[2].DistMax = R;

Background Field = 2;

// --- Generar y guardar --------------------------------------
Mesh 3;
Save "torus_coil.msh";