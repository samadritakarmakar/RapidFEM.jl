 
//+
Point(1) = {0.0, 0.0, 0, 1.0};
//+
Point(2) = {100.0, 0.0, 0, 1.0};
//+
Point(3) = {100.0, 100.0, 0, 1.0};
//+
Point(4) = {0, 100.0, 0, 1.0};
//+
Line(1) = {1, 1};
//+
Line(2) = {1, 1};
//+
Line(3) = {1, 1};
//+
Line(4) = {1, 1};
//+
Line(5) = {1, 2};
//+
Line(6) = {2, 3};
//+
Line(7) = {3, 4};
//+
Line(8) = {4, 1};
//+
Curve Loop(1) = {8, 5, 6, 7};
//+
Plane Surface(1) = {1};
//+
Physical Curve("bottom", 1) = {5};
//+
Physical Curve("top", 2) = {7};
//+
Physical Curve("right", 3) = {6};
//+
Physical Curve("left", 4) = {8};
//+
Physical Surface("vol", 5) = {1};
//+
Transfinite Curve {8, 7, 6, 5} = 10 Using Progression 1;
