SetFactory("OpenCASCADE"); 
//+
Point(1) = {0.0, 0.0, 0, 1.0};
//+
Point(2) = {100.0, 0.0, 0, 1.0};
//+
Point(3) = {100.0, 100.0, 0, 1.0};
//+
Point(4) = {0, 100.0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+

//+

Circle(5) = {50, 50, 0, 25, 0, 2*Pi};
//+
Curve Loop(1) = {5};
//+
Curve Loop(2) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1, 2};
//+
Transfinite Curve {2, 3, 4, 1} = 10 Using Progression 1;
//+
Transfinite Curve {5} = 20 Using Progression 1;
//+
Physical Curve("bottom", 1) = {1};
//+
Physical Curve("top", 2) = {3};
//+
Physical Curve("right", 3) = {2};
//+
Physical Curve("left", 4) = {4};
//+
Physical Surface("vol", 5) = {1};
