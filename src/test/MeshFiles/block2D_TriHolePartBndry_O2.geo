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

//+
Line(1) = {4, 5};
//+
Line(2) = {5, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {4, 1};
//+
Line(5) = {1, 2};
//+
Curve Loop(1) = {1, 2, 3, 4, 5};
//+


Point(6) = {25, 25, 0, 1.0};
//+
Point(7) = {75, 25, 0, 1.0};
//+
Point(8) = {50, 75, 0, 1.0};
//+
Line(9) = {6, 7};
//+
Line(10) = {7, 8};
//+
Line(11) = {8, 6};
//+
Curve Loop(2) = {9, 10, 11};


//+
Plane Surface(1) = {1, 2};
//+
Transfinite Curve {1} = 7 Using Progression 1;
//+
Transfinite Curve {2} = 8 Using Progression 1;
//+
Transfinite Curve {3, 4, 5} = 10 Using Progression 1;

//+
Transfinite Curve {9, 10, 11} = 5 Using Progression 1;

//+
Physical Curve("bottom", 1) = {5};
//+
Physical Curve("top", 2) = {2};
//+
Physical Curve("right", 3) = {3};
//+
Physical Curve("left", 4) = {4};
//+
Physical Surface("vol", 5) = {1};
