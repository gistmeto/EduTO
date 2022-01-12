// Gmsh project created on Mon Jun 21 16:09:35 2021
SetFactory("OpenCASCADE");

Box(1) = {0.1, 0, 0, 0.1, 0.04, 0.04};
Box(2) = {0, 0, 0.04, 0.2, 0.04, 0.002};
Box(3) = {0, 0, 0.042, 0.2, 0.04, 0.04};
Box(4) = {0, 0, 0.042, 0.2, 0.03, 0.03};
BooleanDifference{ Volume{3}; Delete; }{ Volume{4}; Delete; }
Box(5) = {0, 0, -0.005, 0.205, 0.045, 0.092};
BooleanDifference{ Volume{5}; Delete; }{ Volume{1}; Volume{2}; Volume{3}; }

Rectangle(120) = {0.19, 0, 0, 0.01, 0.04, 0};
//+
BooleanFragments{ Volume{1}; Volume{5}; Volume{2}; Volume{3}; Surface{120}; Delete; }{ }
//+
Physical Volume("Design", 144) = {3, 1};
//+
Physical Volume("NonDesign", 145) = {2};
//+
Physical Volume("External", 146) = {5};
//+
Physical Surface("Fixed", 147) = {120};
//+
Physical Surface("Traction", 148) = {137};
//+
Physical Surface("RollerU", 149) = {142, 139};
//+
Physical Surface("RollerV", 150) = {143, 122, 140};
//+
Physical Surface("DesignBoundary", 151) = {137, 22, 24, 135, 27, 130, 123, 121, 124, 120};
