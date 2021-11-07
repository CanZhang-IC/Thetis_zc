res = 20;
t_res = 5;
rad = 27;

x_end = 2000;
y_end = 600;

Point(1) ={0,0,0,res};
Point(2) ={x_end,0,0,res};
Point(3) ={x_end,y_end,0,res};
Point(4) ={0,y_end,0,res};

//circle centre point
Point(5) = {x_end/2,y_end/2,0,t_res};
//left
Point(6) = {x_end/2-rad,y_end/2,0,t_res};
//right
Point(7) = {x_end/2+rad,y_end/2,0,t_res};
//top
Point(8) = {x_end/2,y_end/2+rad,0,t_res};
//bottom
Point(9) = {x_end/2,y_end/2-rad,0,t_res};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Circle(11) = {6,5,7};
Circle(12) = {7,5,6};
Line Loop (10) = {1,2,3,4};
Line Loop (11) = {11,12};
Line Loop (12) = {11,12};

Plane Surface (12) = {10,11};
Plane Surface (13) = {12};

Physical Surface(66) = {12};
Physical Surface(1) = {13};
Physical Line (1) ={4};
Physical Line (2) ={1,3};
Physical Line (3) ={2};
