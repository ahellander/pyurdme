lc = 0.1e-6;
Point(1) = {0, 0, 0, lc};
Point(2) = {0, -.25e-6, 0, lc};
Point(3) = {.25e-6, 0, 0, lc};
Point(4) = {0, 2e-6, 0, lc};
Point(5) = {0, 2.25e-6, 0, lc};
Point(6) = {0.25e-6, 2e-6, 0, lc};

Circle(1) = {5, 4, 6};
Circle(2) = {2, 1, 3};
Line(3) = {6, 3};
Line(4) = {2, 5};
Line Loop(5) = {1, 3, -2, 4};

Plane Surface(6) = {5};

Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Surface{6};
}
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Surface{23};
}
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Surface{40};
}
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Surface{57};
}

Physical Volume(73) = {1, 2, 3, 4};
