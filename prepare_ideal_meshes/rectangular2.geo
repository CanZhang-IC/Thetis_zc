basin_x = 2000;
basin_y = 600;

headland_x_scale = 0.2;
headland_y = 200;

element_size = 4;
element_size_coarse = 40;

Point(1) = {0, 0, 0, element_size_coarse};
Point(2) = {basin_x, 0, 0, element_size_coarse};
Point(3) = {basin_x, basin_y, 0, element_size_coarse};
Point(4) = {0, basin_y, 0, element_size_coarse};


Line(100) = {3, 4};
Line(101) = {4, 1};
Line(102) = {1, 2};
Line(103) = {2, 3};
Line Loop(104) = {100, 101, 102, 103};

// Generate site nodes
Point(1000) = {750, 150, 0, element_size};
Extrude{500, 0, 0} { Point{1000}; Layers{500/element_size}; }
Extrude{0, 300, 0} { Line{105}; Layers{300/element_size}; }
Line Loop(110) = {106, -108, -105, 107};
Plane Surface(111) = {104, 110};
Physical Line(1) = {101};
Physical Line(2) = {103};
Physical Line(3) = {100, 102};
Physical Surface(1) = {111};
Physical Surface(2) = {109};

