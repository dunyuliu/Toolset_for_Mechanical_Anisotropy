// Gmsh project created on Thu Apr 28 12:03:24 2022
SetFactory("OpenCASCADE");
dx_schist = 0.1;
dx_other_rock = 1;
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 10, 0, 1.0};
//+
Point(3) = {10, 10, 0, 1.0};
//+
Point(4) = {10, 0, 0, 1.0};
//+
Point(5) = {3, 0, -3, 1.0};
//+
Point(6) = {3, 10, -3, 1.0};
//+
Point(7) = {10, 0, -3, 1.0};
//+
Point(8) = {10, 10, -3, 1.0};
//+
Point(9) = {3.5, 3, 0, dx_schist};
//+
Point(10) = {7.8, 5.5, 0, dx_schist};
//+
Point(11) = {7.3, 6.4, 0, dx_schist};
//+
Point(12) = {3, 3.9, 0, dx_schist};
//+
Point(13) = {3, 3.9, -1, dx_schist};
//+
Point(14) = {7.3, 6.4, -1, dx_schist};
//+
Line(1) = {1, 5};
//+
Line(2) = {5, 6};
//+
Line(3) = {6, 2};
//+
Line(4) = {2, 1};
//+
Line(5) = {2, 3};
//+
Line(6) = {3, 8};
//+
Line(7) = {8, 6};
//+
Line(8) = {5, 7};
//+
Line(9) = {7, 4};
//+
Line(10) = {4, 3};
//+
Line(11) = {8, 7};
//+
Line(12) = {12, 9};
//+
Line(13) = {9, 13};
//+
Line(14) = {13, 12};
//+
Line(15) = {12, 11};
//+
Line(16) = {11, 10};
//+
Line(17) = {10, 9};
//+
Line(18) = {13, 14};
//+
Line(19) = {14, 11};
//+
Line(20) = {10, 14};
//+
Line(21) = {1, 4};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, 7, 3};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, 6, 11, 9};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {21, -9, -8, -1};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {2, -7, 11, -8};
//+
Plane Surface(5) = {5};
//+
Line(22) = {2, 12};
//+
Line(23) = {9, 1};
//+
Line(24) = {10, 4};
//+
Line(25) = {11, 3};
//+
Curve Loop(6) = {22, 12, 23, -4};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {15, 25, -5, 22};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {16, 24, 10, -25};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {17, 23, 21, -24};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {18, -20, 17, 13};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {16, 20, 19};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {15, 16, 17, -12};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {18, 19, -15, -14};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {12, 13, 14};
//+
Plane Surface(14) = {14};
//+
Surface Loop(1) = {12, 13, 10, 11, 14};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {7, 8, 9, 6, 1, 4, 3, 2, 5, 13, 10, 11, 14};
//+
Volume(2) = {2};
//+
Physical Surface("south", 26) = {4};
//+
Physical Surface("north", 27) = {2};
//+
Physical Surface("west", 28) = {1};
//+
Physical Surface("east", 29) = {3};
//+
Physical Surface("bottom", 30) = {5};
//+
Physical Volume("schist", 31) = {1};
//+
Physical Volume("other", 32) = {2};
