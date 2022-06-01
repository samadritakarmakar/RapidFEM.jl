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
Point(5) = {70, 100.0, 0, 1.0};
//+
Circle(5) = {50, 50, 0, 25, 0, 2*Pi};
//+

//+
Line(6) = {4, 5};
//+
Line(7) = {5, 3};
//+
Line(8) = {3, 2};
//+
Line(9) = {4, 1};
//+
Line(10) = {1, 2};
//+
Curve Loop(1) = {6, 7, 8, -10, -9};
//+
Curve Loop(2) = {5};
//+
Plane Surface(1) = {1, 2};
//+
Transfinite Curve {6} = 7 Using Progression 1;
//+
Transfinite Curve {7} = 8 Using Progression 1;
//+
Transfinite Curve {10, 9} = 10 Using Progression 1;
//+
Transfinite Curve {8} = 15 Using Progression 1;
//+
Transfinite Curve {5} = 20 Using Progression 1;

//+
Physical Curve("bottom", 1) = {10};
//+
Physical Curve("top", 2) = {7};
//+
Physical Curve("right", 3) = {8};
//+
Physical Curve("left", 4) = {9};
//+
Physical Surface("vol", 5) = {1};
