// ============================================================
// Bobina de Helmholtz TRIAXIAL CUADRADA — 3 pares ortogonales
// Bobinas filamentarias rectas (4 lados por cuadrado)
// ============================================================
SetFactory("OpenCASCADE");

// --- Parámetro externo (barrible desde Julia) -----------------
// Uso: onelab.setNumber("Parameters/L", [L]) antes de abrir el .geo
L = DefineNumber[ 1.00, Name "Parameters/L" ];   // lado BASE (bobina Z) [m]

// --- Factores de anidado por eje (evitan intersección física) --
fx = 1.00;
fy = 1.08;
fz = 1.16;

Lx = L * fx;
Ly = L * fy;
Lz = L * fz;

// --- Separación Helmholtz ÓPTIMA para bobina CUADRADA -----------
// (a diferencia de la circular, aquí d = 0.5445*L, no L/2)
kopt = 0.5445;
dx = kopt * Lx;
dy = kopt * Ly;
dz = kopt * Lz;

// --- Dominio de aire (esfera) ------------------------------------
// "Alcance" característico = medio lado de la bobina más grande
Rext  = Lz / 2;
R_dom = 3.0 * Rext;

// --- Ancho físico del hilo/fuente (debe coincidir con σ_frac en Julia) --
// Proporcional a L (no un valor absoluto fijo): mantiene la resolución
// relativa constante si luego barres distintos L. El valor L/40 (25mm
// para L=1m) se escogió por robustez de malla, no por un gauge de
// alambre real: valores más finos (p.ej. L/100) generan una malla con
// calidad promedio muy pobre (<0.05) porque el dominio de aire es
// ~1.7 m de radio y el salto de tamaño resulta demasiado brusco para
// que Delaunay lo gradúe bien (ver dist_min/dist_max más abajo).
sigma_src = L / 40;

lc_coil_L     = L / 20;
lc_coil_sigma = sigma_src;
If (lc_coil_L < lc_coil_sigma)
  lc_coil = lc_coil_L;
Else
  lc_coil = lc_coil_sigma;
EndIf
lc_far = R_dom / 5;

Sphere(1) = {0, 0, 0, R_dom};

// ============================================================
// --- Puntos de las 6 bobinas cuadradas --------------------------
// Orden de los 4 vértices: sentido CCW visto desde el lado + del
// eje de esa bobina (mismo orden en las 2 bobinas del par, para
// que la corriente circule igual y el campo se sume en el centro)
// ============================================================

// Par Z: cuadrado en plano XY, en z=+dz y z=-dz
Point(501) = {-Lz/2, -Lz/2,  dz};
Point(502) = { Lz/2, -Lz/2,  dz};
Point(503) = { Lz/2,  Lz/2,  dz};
Point(504) = {-Lz/2,  Lz/2,  dz};

Point(511) = {-Lz/2, -Lz/2, -dz};
Point(512) = { Lz/2, -Lz/2, -dz};
Point(513) = { Lz/2,  Lz/2, -dz};
Point(514) = {-Lz/2,  Lz/2, -dz};

// Par X: cuadrado en plano YZ, en x=+dx y x=-dx
Point(521) = { dx, -Lx/2, -Lx/2};
Point(522) = { dx,  Lx/2, -Lx/2};
Point(523) = { dx,  Lx/2,  Lx/2};
Point(524) = { dx, -Lx/2,  Lx/2};

Point(531) = {-dx, -Lx/2, -Lx/2};
Point(532) = {-dx,  Lx/2, -Lx/2};
Point(533) = {-dx,  Lx/2,  Lx/2};
Point(534) = {-dx, -Lx/2,  Lx/2};

// Par Y: cuadrado en plano XZ, en y=+dy y y=-dy
Point(541) = {-Ly/2, dy, -Ly/2};
Point(542) = {-Ly/2, dy,  Ly/2};
Point(543) = { Ly/2, dy,  Ly/2};
Point(544) = { Ly/2, dy, -Ly/2};

Point(551) = {-Ly/2, -dy, -Ly/2};
Point(552) = {-Ly/2, -dy,  Ly/2};
Point(553) = { Ly/2, -dy,  Ly/2};
Point(554) = { Ly/2, -dy, -Ly/2};

// --- Líneas (4 lados por cuadrado, cerrando el lazo) ------------
Line(201) = {501,502}; Line(202) = {502,503}; Line(203) = {503,504}; Line(204) = {504,501};
Line(211) = {511,512}; Line(212) = {512,513}; Line(213) = {513,514}; Line(214) = {514,511};

Line(221) = {521,522}; Line(222) = {522,523}; Line(223) = {523,524}; Line(224) = {524,521};
Line(231) = {531,532}; Line(232) = {532,533}; Line(233) = {533,534}; Line(234) = {534,531};

Line(241) = {541,542}; Line(242) = {542,543}; Line(243) = {543,544}; Line(244) = {544,541};
Line(251) = {551,552}; Line(252) = {552,553}; Line(253) = {553,554}; Line(254) = {554,551};

// --- Embeber las 24 líneas en el volumen de aire -----------------
Curve{201:204, 211:214, 221:224, 231:234, 241:244, 251:254} In Volume{1};

// --- Grupos físicos ------------------------------------------------
Physical Volume("air")       = {1};
Physical Curve("coilZ1")     = {201:204};
Physical Curve("coilZ2")     = {211:214};
Physical Curve("coilX1")     = {221:224};
Physical Curve("coilX2")     = {231:234};
Physical Curve("coilY1")     = {241:244};
Physical Curve("coilY2")     = {251:254};
Physical Surface("boundary") = {1};

// --- Delaunay clásico (evita el bucle costoso de HXT) ------------
Mesh.Algorithm3D = 1;

// --- Refinamiento cerca de las 24 líneas --------------------------
Field[1] = Distance;
Field[1].CurvesList = {201:204, 211:214, 221:224, 231:234, 241:244, 251:254};
Field[1].Sampling   = 150;

// Zona de transición atada a sigma_src (que ya escala con L): da
// una banda fina alrededor del hilo y deja que el resto del campo
// (Threshold hasta lc_far) gradúe el tamaño en el resto del dominio.
dist_min = 1.5 * sigma_src;
dist_max = 5   * sigma_src;

Field[2] = Threshold;
Field[2].InField  = 1;
Field[2].SizeMin  = lc_coil;
Field[2].SizeMax  = lc_far;
Field[2].DistMin  = dist_min;
Field[2].DistMax  = dist_max;

// Refinar la región central (zona de medición del campo uniforme)
Field[3] = Ball;
Field[3].Radius    = Rext / 2;
Field[3].Thickness = Rext / 4;
Field[3].VIn  = Rext / 15;
Field[3].VOut = lc_far;
Field[3].XCenter = 0; Field[3].YCenter = 0; Field[3].ZCenter = 0;

Field[10] = Min;
Field[10].FieldsList = {2, 3};
Background Field = 10;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints         = 0;
Mesh.MeshSizeFromCurvature      = 0;
