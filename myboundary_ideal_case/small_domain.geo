element_size = 150;
element_size_coarse = 400;

Point(1) = {5000, 1000, 0, element_size_coarse};
Point(2) = {15000, 1000, 0, element_size_coarse};
Point(3) = {15000, 3000, 0, element_size_coarse};
Point(4) = {5000, 3000, 0, element_size_coarse};



Line(100) = {3, 4};
Line(101) = {4, 1};
Line(102) = {1, 2};
Line(103) = {2, 3};
Line Loop(104) = {100, 101, 102, 103};


Plane Surface(111) = {104};

Physical Line(1) = {101};
Physical Line(2) = {103};
Physical Line(3) = {100, 102};
Physical Surface(1) = {111};

