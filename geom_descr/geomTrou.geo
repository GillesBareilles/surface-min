h = 0.1;

Point(1) = {0, 0, 0, h};
Point(2) = {0.5, 0, 0, h};
Point(3) = {2, 0, 0, h};
Point(4) = {-0.5, 0, 0, h};
Point(5) = {-2, 0, 0, h};

Circle(1) = {2, 1, 4};
Circle(2) = {3, 1, 5};
Circle(3) = {4, 1, 2};
Circle(4) = {5, 1, 3};
Line Loop(5) = {2, 4};
Line Loop(6) = {1, 3};
Plane Surface(7) = {5, 6};
