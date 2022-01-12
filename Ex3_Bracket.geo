Merge "Ex3_Bracket.STEP";
//+
SetFactory("OpenCASCADE");
Box(5) = {-12.5, -20, 0, 132.5, 90, 25};
//+
BooleanDifference{Volume{5}; Delete; }{Volume{4}; Volume{3}; Volume{2}; Volume{1}; }
//+
BooleanFragments{Volume{5};Volume{4}; Volume{3}; Volume{2}; Volume{1}; Delete; }{ }
//+
Physical Volume("Design", 141) = {3};
//+
Physical Volume("NonDesign", 142) = {1, 2, 4};
//+
Physical Volume("External", 143) = {5};
//+
Physical Surface("Fixed", 144) = {4, 10};
//+
Physical Surface("Traction", 145) = {43};
//+
Physical Surface("RollerW", 146) = {66, 71};
