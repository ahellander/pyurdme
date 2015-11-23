// Gmsh project created on Wed Aug  7 10:53:48 2013
Point(1) = {0, 0.5e-6, 0, 0.5e-6};
Point(2) = {0.5e-6, 0.5e-6, 0, 0.5e-6};
Point(3) = {0, 0, 0, 0.5e-6};
Point(4) = {0.5e-6, 0, 0, 0.5e-6};
Line(1) = {1, 2};
Line(2) = {1, 3};
Line(3) = {3, 4};
Line(4) = {4, 2};
Line Loop(5) = {4, -1, 2, 3};
Plane Surface(6) = {5};
Physical Surface(7) = {6};
Physical Line(8) = {2, 1, 4, 3};
Coherence;
