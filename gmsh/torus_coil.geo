// ============================================================
// Bobina de Helmholtz — bobinas filamentarias (curvas 1D)
// Malla volumétrica de aire + curvas embebidas como fuente
// ============================================================
SetFactory("OpenCASCADE");

// --- Parámetros ---------------------------------------------
R     = 0.60;      // radio [m] (60 cm)
sep   = 0.60;      // separación entre planos [m] (60 cm): par de Helmholtz
R_dom = 3.0 * R;   // radio del dominio de aire

lc_coil = R / 40;  // tamaño de malla sobre las bobinas (fino)
lc_far  = R_dom/5; // tamaño de malla en la frontera (grueso)

// --- Dominio de aire (esfera) -------------------------------
Sphere(1) = {0, 0, 0, R_dom};

// --- Bobinas como círculos (curvas) -------------------------
// Planos z = ±sep/2. La distancia entre bobinas es sep = R.
Circle(100) = {0, 0,  sep/2, R, 0, 2*Pi};
Circle(101) = {0, 0, -sep/2, R, 0, 2*Pi};

// --- Control de tamaño de malla -----------------------------
Mesh.Algorithm3D = 10;   // HXT: rápido y robusto para volúmenes grandes

// Refinar cerca de las curvas de las bobinas
Field[1] = Distance;
Field[1].CurvesList = {100, 101};
Field[1].Sampling   = 200;

Field[2] = Threshold;
Field[2].InField  = 1;
Field[2].SizeMin  = lc_coil;
Field[2].SizeMax  = lc_far;
Field[2].DistMin  = R / 10;   // zona fina alrededor del hilo
Field[2].DistMax  = R;        // transición suave hacia el exterior

// Refinar también la región central (donde medimos el campo uniforme)
Field[3] = Ball;
Field[3].Radius   = R / 2;
Field[3].Thickness = R / 4;
Field[3].VIn  = R / 25;       // fino en el centro de Helmholtz
Field[3].VOut = lc_far;
Field[3].XCenter = 0; Field[3].YCenter = 0; Field[3].ZCenter = 0;

Field[10] = Min;
Field[10].FieldsList = {2, 3};
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
