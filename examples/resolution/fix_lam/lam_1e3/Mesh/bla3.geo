cl1 = 1./8.;
Point(1) = {-1, 0, 0, cl1};
Point(2) = {1, 0, 0, cl1};
Point(3) = {-1, -1, 0, cl1};
Point(4) = {1, -1, 0, cl1};
Point(5) = {0, -0, 0, cl1};
Point(6) = {0, -0.2, 0, cl1};
Point(7) = {0, -0.4, 0, cl1};
Point(8) = {0, -0.25, 0, cl1};
Point(9) = {0, -0.3, 0, cl1};
Point(10) = {0, -0.35, 0, cl1};
Line(1) = {3, 1};
Line(2) = {1, 5};
Line(3) = {5, 2};
Line(4) = {2, 4};
Line(5) = {4, 3};
Line(6) = {6, 7};
Line Loop(7) = {1, 2, 3, 4,5};
Plane Surface(8) = {7};
Point{10, 9, 8, 7, 6} In Surface{8};
