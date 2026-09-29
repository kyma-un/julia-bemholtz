// ============================================================
// Bobina de Helmholtz — bobinas filamentarias (curvas 1D)
// Malla volumétrica de aire + curvas embebidas como fuente
// ============================================================
SetFactory("OpenCASCADE");

// --- Parámetros ---------------------------------------------
R     = 0.60;      // radio [m] (60 cm)
sep   = 0.60;      // separación entre planos [m] (60 cm): par de Helmholtz
R_dom = 3.0 * R;   // radio del dominio de aire

lc_coil = R / 36;   // tamaño de malla sobre las bobinas
lc_in   = 0.022;     // tamaño en el centro
lc_far  = 0.45;      // tamaño en la frontera de la esfera

// --- Dominio de aire (esfera) -------------------------------
Sphere(1) = {0, 0, 0, R_dom};

// --- Bobinas como círculos (curvas) -------------------------
// Planos z = ±sep/2. La distancia entre bobinas es sep = R.
Circle(100) = {0, 0,  sep/2, R, 0, 2*Pi};
Circle(101) = {0, 0, -sep/2, R, 0, 2*Pi};

// --- Control de tamaño de malla -----------------------------
Mesh.Algorithm3D = 1;    // Delaunay: tetraedros más parejos

// h = h(r): el mismo tamaño a un radio dado, en cualquier dirección.
// Crece con r*r para no dejar el interior a 3 cm hasta la frontera.
Field[6] = MathEval;
Field[6].F = "0.022 + 0.428*(Sqrt(x*x+y*y+z*z)/1.8)*(Sqrt(x*x+y*y+z*z)/1.8)";

// Refinar cerca de las curvas de las bobinas
Field[1] = Distance;
Field[1].CurvesList = {100, 101};
Field[1].Sampling   = 200;

Field[2] = Threshold;
Field[2].InField  = 1;
Field[2].SizeMin  = lc_coil;
Field[2].SizeMax  = lc_far;
Field[2].DistMin  = 0.035;    // zona fina alrededor del hilo
Field[2].DistMax  = 0.14;     // el resto lo da el campo radial

Field[10] = Min;
Field[10].FieldsList = {2, 6};
Background Field = 10;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// Si otro .geo incluye este archivo, define AS_INCLUDE = 1 antes del Include
// y se queda solo con la geometría. Abierto solo, mallar como siempre.
If (Exists(AS_INCLUDE) == 0)
  // Embeber las curvas en el aire para que la malla las respete
  Curve{100, 101} In Volume{1};

  Physical Volume("air")       = {1};
  Physical Curve("coil1")      = {100};
  Physical Curve("coil2")      = {101};
  Physical Surface("boundary") = {1};

  Mesh 3;
  Save "torus_coil.msh";
EndIf
