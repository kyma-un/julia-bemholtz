// ============================================================
// Magnetorquers — dos solenoides ortogonales
// Solo estas bobinas. El aire y Helmholtz están en otros .geo
// ============================================================
SetFactory("OpenCASCADE");

// --- Parámetros ---------------------------------------------
If (Exists(mt_diameter) == 0)
  mt_diameter = 0.010;   // diámetro [m]  (10 mm)
EndIf
If (Exists(mt_length) == 0)
  mt_length = 0.140;     // largo [m]      (140 mm)
EndIf

mt_gap    = 0.003;                 // holgura para que no se toquen [m]
mt_radius = mt_diameter / 2;

// Las bobinas de Helmholtz están en planos horizontales (normal +Z).
// Estos solenoides viven en el plano XZ (y = 0), perpendicular a esos planos.
// Una a lo largo de X y la otra a lo largo de Z, centradas en el origen.
z_join = mt_length / 2 + mt_radius + mt_gap;
z_mid  = z_join / 2;

// Bobina en X
Cylinder(10) = {-mt_length/2, 0, z_mid, mt_length, 0, 0, mt_radius};
// Bobina en Z, ortogonal
Cylinder(11) = {0, 0, -mt_length/2 - z_mid, 0, 0, mt_length, mt_radius};

If (Exists(AS_INCLUDE) == 0)
  Physical Volume("magnetorquer_x") = {10};
  Physical Volume("magnetorquer_z") = {11};

  Mesh.CharacteristicLengthMin = mt_diameter / 12;
  Mesh.CharacteristicLengthMax = mt_diameter / 6;
  Mesh 3;
  Save "magnetorquers.msh";
EndIf
