// ============================================================
// Bobina de Helmholtz TRIAXIAL — 3 pares ortogonales (X, Y, Z)
// Bobinas filamentarias (curvas 1D) + malla volumétrica de aire
// ============================================================
SetFactory("OpenCASCADE");

// --- Parámetro externo (barrible desde Julia/CLI) -----------
// Uso: gmsh triaxial_coil.geo -setnumber R 0.08 -3 -o out.msh
R = DefineNumber[ 0.100, Name "Parameters/R" ];   // radio BASE [m]

// --- Factores de anidado por eje -----------------------------
// Evitan que las 3 bobinas se intersecten físicamente (nesting).
// Ajusta si tu diseño real usa otra separación entre radios.
fx = 1.00;
fy = 1.08;
fz = 1.16;

Rx = R * fx;
Ry = R * fy;
Rz = R * fz;

// Separación Helmholtz de cada par: d = R_eje/2
dx = Rx / 2;
dy = Ry / 2;
dz = Rz / 2;

// --- Dominio de aire (esfera) ---------------------------------
// Debe contener la bobina más grande (Rz) con margen.
R_dom = 3.0 * Rz;

lc_coil = R / 20;   // tamaño de malla sobre las bobinas (fino, ajustado para no agotar RAM)
lc_far  = R_dom/5;  // tamaño de malla en la frontera (grueso)

Sphere(1) = {0, 0, 0, R_dom};

// ============================================================
// --- Par Z: círculos en plano XY, separados en z -------------
// (idéntico a tu geo uniaxial original)
// ============================================================
Circle(100) = {0, 0,  dz, Rz, 0, 2*Pi};   // z+
Circle(101) = {0, 0, -dz, Rz, 0, 2*Pi};   // z-

// ============================================================
// --- Par X: círculos en plano YZ, separados en x -------------
// Se crean en plano XY (normal z) y se rotan 90° sobre eje Y
// para que su normal quede alineada con X.
// ============================================================
Circle(102) = {0, 0, 0, Rx, 0, 2*Pi};
Rotate{ {0,1,0}, {0,0,0}, Pi/2 } { Curve{102}; }
Translate{ dx, 0, 0 } { Curve{102}; }     // x+

Circle(103) = {0, 0, 0, Rx, 0, 2*Pi};
Rotate{ {0,1,0}, {0,0,0}, Pi/2 } { Curve{103}; }
Translate{ -dx, 0, 0 } { Curve{103}; }    // x-

// ============================================================
// --- Par Y: círculos en plano XZ, separados en y -------------
// Se rotan 90° sobre eje X para que su normal quede en Y.
// (misma rotación en ambos -> mismo sentido de corriente)
// ============================================================
Circle(104) = {0, 0, 0, Ry, 0, 2*Pi};
Rotate{ {1,0,0}, {0,0,0}, Pi/2 } { Curve{104}; }
Translate{ 0, dy, 0 } { Curve{104}; }     // y+

Circle(105) = {0, 0, 0, Ry, 0, 2*Pi};
Rotate{ {1,0,0}, {0,0,0}, Pi/2 } { Curve{105}; }
Translate{ 0, -dy, 0 } { Curve{105}; }    // y-

// --- Embeber las 6 curvas en el volumen de aire --------------
Curve{100, 101, 102, 103, 104, 105} In Volume{1};

// --- Grupos físicos -------------------------------------------
Physical Volume("air")       = {1};
Physical Curve("coilZ1")     = {100};
Physical Curve("coilZ2")     = {101};
Physical Curve("coilX1")     = {102};
Physical Curve("coilX2")     = {103};
Physical Curve("coilY1")     = {104};
Physical Curve("coilY2")     = {105};
Physical Surface("boundary") = {1};   // superficie exterior de la esfera

// --- Control de tamaño de malla ---------------------------------
Mesh.Algorithm3D = 10;   // HXT: rápido y robusto para volúmenes grandes

// Refinar cerca de las 6 curvas de las bobinas
Field[1] = Distance;
Field[1].CurvesList = {100, 101, 102, 103, 104, 105};
Field[1].Sampling   = 150;

Field[2] = Threshold;
Field[2].InField  = 1;
Field[2].SizeMin  = lc_coil;
Field[2].SizeMax  = lc_far;
Field[2].DistMin  = R / 10;   // zona fina alrededor del hilo
Field[2].DistMax  = R;        // transición suave hacia el exterior

// Refinar también la región central (donde medimos el campo uniforme)
Field[3] = Ball;
Field[3].Radius    = R / 2;
Field[3].Thickness = R / 4;
Field[3].VIn  = R / 15;       // fino en el centro de Helmholtz (ajustado)
Field[3].VOut = lc_far;
Field[3].XCenter = 0; Field[3].YCenter = 0; Field[3].ZCenter = 0;

Field[10] = Min;
Field[10].FieldsList = {2, 3};
Background Field = 10;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints         = 0;
Mesh.MeshSizeFromCurvature      = 0;

// NOTA: no se llama "Mesh 3;" ni "Save" aquí.
// El control de mallado/escritura lo hace Julia (Gmsh.jl),
// para poder barrer R en un loop sin regenerar el .geo.
