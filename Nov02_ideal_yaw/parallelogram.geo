basin_x = 2000;
basin_y = 600;

headland_x_scale = 0.2;
headland_y = 200;

site_x = 160;
site_y = 100;
site_x_start = basin_x/2-site_x/2;
site_x_end = basin_x/2+site_x/2;

site_y_start = basin_y/2 - 50;
site_y_end = site_y_start+site_y;

element_size = 4;
element_size_coarse = 20;

Point(1) = {600, 0, 0, element_size_coarse};
Point(2) = {1000, 0, 0, element_size_coarse};
Point(3) = {1000, basin_y, 0, element_size_coarse};
Point(4) = {1400, basin_y, 0, element_size_coarse};

Line(100) = {4, 3};
Line(101) = {3, 1};
Line(102) = {1, 2};
Line(103) = {2, 4};
Line Loop(104) = {100, 101, 102, 103};

// Generate site nodes
Point(1000) = {site_x_start, site_y_start, 0, element_size};
Extrude{site_x, 0, 0} { Point{1000}; Layers{site_x/element_size}; }
Extrude{0, site_y, 0} { Line{105}; Layers{site_y/element_size}; }
Line Loop(110) = {106, -108, -105, 107};
Plane Surface(111) = {104, 110};
Physical Line(1) = {100};
Physical Line(2) = {102};
Physical Line(3) = {101, 103};
Physical Surface(1) = {111};
Physical Surface(2) = {109};
