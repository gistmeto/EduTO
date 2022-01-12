// Gmsh project created on Mon Oct 04 16:36:22 2021
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, 0, 0.3, 0, 0, 0.1, Pi/2};
Rotate {{1, 0, 0}, {0, 0, 0}, -Pi/2} {Volume{1};}

Point(11) = {0.01, 0, 0.1, 1.0};
Point(23) = {0.01, 0.1*Sin(Pi/180*10), 0.1*Cos(Pi/180*10), 1.0};
Point(24) = {0, 0.1*Sin(Pi/180*10), 0.1*Cos(Pi/180*10), 1.0};
Point(25) = {0.01, 0, 0, 1.0};

Line(10) = {3, 11};
Line(11) = {24, 23};
Circle(12) = {3, 6, 24};
Circle(13) = {11, 25, 23};

Rectangle(6) = {0.29, 0.08, 0, 0.01, 0.02, 0};
Curve Loop(7) = {11, -13, -10, 12};
Surface(7) = {7};

Cylinder(100) = {0, 0, 0, 0.31, 0, 0, 0.11, Pi/2};
Rotate {{1, 0, 0}, {0, 0, 0}, -Pi/2} {Volume{100};}
Box(101) = {0, 0, -0.01, 0.31, 0.06, 0.01};

Dilate {{0, 0, 0}, {1, 0.5, 1}} {Volume{1}; Surface{6}; Surface{7};}
Dilate {{0, 0, 0}, {1, 0.5454, 1}} {Volume{100};}

BooleanDifference{ Volume{100}; Delete; }{ Volume{1}; }
BooleanUnion{ Volume{100}; Delete; }{ Volume{101}; Delete; }
BooleanFragments{ Volume{2}; Volume{1}; Surface{7}; Surface{6}; Delete; }{ }

Physical Volume("Design", 86) = {1};
Physical Volume("External", 87) = {2};
Physical Surface("RollerU", 88) = {18};
Physical Surface("RollerV", 89) = {19};
Physical Surface("RollerW", 92) = {6};
Physical Surface("Traction", 91) = {7};

Physical Surface("DesignBoundary", 93) = {7, 13, 14, 17, 6};
