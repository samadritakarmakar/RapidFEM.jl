 
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
Point(5) = {25, 25, 0, 1.0};
//+
Point(6) = {75, 25, 0, 1.0};
//+
Point(7) = {50, 75, 0, 1.0};

//+
Line(9) = {5, 6};
//+
Line(10) = {6, 7};
//+
Line(11) = {7, 5};
//+
Curve Loop(2) = {11, 9, 10};
//+
Plane Surface(1) = {1, 2};
//+
Transfinite Curve {5, 6, 7, 8} = 10 Using Progression 1;
//+
Transfinite Curve {9, 10, 11} = 5 Using Progression 1;
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
