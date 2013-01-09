lc = 0.1;
Point(1) = {0, 0, 0, lc};
Point(2) = {0, -.25, 0, lc};
Point(3) = {.25, 0, 0, lc};
Point(4) = {0, 2, 0, lc};
Point(5) = {0, 2.25, 0, lc};
Point(6) = {0.25, 2, 0, lc};

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

Physical Surface(73) = {31, 48, 65, 14};
Physical Surface(74) = {35, 18, 69, 52};
Physical Surface(75) = {38, 55, 72, 21};
Physical Volume(76)  = {2, 1, 4, 3};
