//+
SetFactory("OpenCASCADE");

// x y z r R giro 
Torus(1) = {0, 0, -0.35, 0.35, 0.01, 2*Pi};
Torus(2) = {0, 0, 0.35, 0.35, 0.01, 2*Pi};

// x y z r1 r2
Box(3) = {-0.5, -0.5, -0.5, 1, 1, 1};