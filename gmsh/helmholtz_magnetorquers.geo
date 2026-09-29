// ============================================================
// Helmholtz + magnetorquers
// Une gmsh/torus_coil.geo y gmsh/magnetorquers.geo
//
// Diámetro y largo de los magnetorquers:
// ============================================================
mt_diameter = 0.010;   // [m]  10 mm
mt_length   = 0.140;   // [m]  140 mm

AS_INCLUDE = 1;
Include "torus_coil.geo";
Include "magnetorquers.geo";

// Mete los dos cilindros dentro de la esfera de aire
BooleanFragments{ Volume{1}; Delete; }{ Volume{10, 11}; Delete; }

eps = 1e-4;
mt_x() = Volume In BoundingBox{
  -mt_length/2 - eps, -mt_radius - eps, z_mid - mt_radius - eps,
   mt_length/2 + eps,  mt_radius + eps, z_mid + mt_radius + eps};
mt_z() = Volume In BoundingBox{
  -mt_radius - eps, -mt_radius - eps, -z_mid - mt_length/2 - eps,
   mt_radius + eps,  mt_radius + eps, -z_mid + mt_length/2 + eps};

air() = Volume{:};
air() -= mt_x();
air() -= mt_z();

If (#air() != 1 || #mt_x() != 1 || #mt_z() != 1)
  Error("No se pudieron separar el aire y los dos magnetorquers.");
EndIf

Curve{100, 101} In Volume{air(0)};

outer() = Boundary{ Volume{air()}; };
mt_skin() = Boundary{ Volume{mt_x(), mt_z()}; };
outer() -= mt_skin();

Physical Volume("air")             = {air()};
Physical Volume("magnetorquer_x")  = {mt_x()};
Physical Volume("magnetorquer_z")  = {mt_z()};
Physical Curve("coil1")            = {100};
Physical Curve("coil2")            = {101};
Physical Surface("boundary")       = {outer()};

// Refino sobre los magnetorquers. El campo 10 llega de torus_coil.geo
Field[4] = Distance;
Field[4].SurfacesList = {mt_skin()};
Field[4].Sampling = 80;

Field[5] = Threshold;
Field[5].InField  = 4;
Field[5].SizeMin  = mt_diameter / 3;
Field[5].SizeMax  = lc_far;   // lejos de las varillas no aplasta el campo radial
Field[5].DistMin  = mt_diameter;
Field[5].DistMax  = 0.05;

Field[11] = Min;
Field[11].FieldsList = {10, 5};
Background Field = 11;

Printf("Magnetorquers: D = %g m, L = %g m", mt_diameter, mt_length);

Mesh 3;
Save "helmholtz_magnetorquers.msh";
