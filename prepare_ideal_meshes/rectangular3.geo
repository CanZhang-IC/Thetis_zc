basin_x = 1000;
basin_y = 100;

site_x = 800;
site_y = 50;

site_x_start = basin_x/2 - site_x/2;
site_y_start = basin_y/2 - site_y/2;

element_size = 2;
element_size_coarse = 10;

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

Point(200) = {0,550,0, element_size_coarse};
Point(201) = {1000,550,0, element_size_coarse};

Line(200) = {4,200};
Line(201) = {200,201};
Line(202) = {201,3};
Line Loop(105) = {200,201,202,100};

Point(300) = {0,-450,0, element_size_coarse};
Point(301) = {1000,-450,0, element_size_coarse};

Line(300) = {1,300};
Line(301) = {300,301};
Line(302) = {301,2};
Line Loop(106) = {300,301,302,-102};

Plane Surface(111) = {104,110};
Plane Surface(112) = {105};
Plane Surface(113) = {106};

Physical Line(1) = {101,200,300};
Physical Line(2) = {103,202,302};
Physical Line(3) = {201, 301};


// Transfinite Line{200,202,300,302,101,103} = 10 Using Progression 1;
// Transfinite Line{201,100,102,301} = 100 Using Progression 1;
// Transfinite Surface(111) = {100, 101, 102, 103};
// Transfinite Surface(105) = {200,201,202,100};
// Transfinite Surface(106) = {300,301,302,-102};

Physical Surface(1) = {111};
Physical Surface(2) = {109};
Physical Surface(3) = {112};
Physical Surface(4) = {113};