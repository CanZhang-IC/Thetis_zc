basin_x = 2000;
basin_y = 600;

site_x = 400;
site_y = 200;

site_x_start = 1400;
site_y_start = basin_y/2 - site_y/2;

element_size = 5;
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
Point(1000) = {site_x_start, site_y_start, 0, element_size};
Extrude{site_x, 0, 0} { Point{1000}; Layers{site_x/element_size}; }
Extrude{0, site_y, 0} { Line{105}; Layers{site_y/element_size}; }
Line Loop(110) = {106, -108, -105, 107};

point_x = 450;
point_y = 250;
area_circle = 40;
//effected area
Point(21) = {point_x - area_circle, point_y - area_circle, 0, element_size_coarse};
Point(22) = {point_x + area_circle, point_y - area_circle, 0, element_size_coarse};
Point(23) = {point_x + area_circle, point_y + area_circle, 0, element_size_coarse};
Point(24) = {point_x - area_circle, point_y + area_circle, 0, element_size_coarse};


Line(2100) = {23, 24};
Line(2101) = {24, 21};
Line(2102) = {21, 22};
Line(2103) = {22, 23};
Line Loop(2104) = {2100, 2101, 2102, 2103};



Plane Surface(111) = {104, 110};
Physical Line(1) = {101};
Physical Line(2) = {103};
Physical Line(3) = {100, 102};
Physical Surface(1) = {111};
Physical Surface(2) = {109};
