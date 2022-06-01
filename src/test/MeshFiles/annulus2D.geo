 
//+
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 3, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 1.5, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Curve Loop(2) = {2};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("outerRadius", 1) = {1};
//+
Physical Curve("innerRadius", 2) = {2};
//+
Transfinite Curve {1} = 30 Using Progression 1;
//+
Transfinite Curve {2} = 20 Using Progression 1;
//+
Physical Surface("vol", 3) = {1};
