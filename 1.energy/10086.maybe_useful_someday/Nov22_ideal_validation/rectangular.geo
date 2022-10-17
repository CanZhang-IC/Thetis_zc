res = 0.2;
t_res = 0.05;
rad = 0.27;

Point(1) ={0,0,0,res};
Point(2) ={50,0,0,res};
Point(3) ={50,1.2,0,res};
Point(4) ={0,1.2,0,res};

//circle centre point
Point(5) = {25,0.6,0,t_res};
//left
Point(6) = {25-rad,0.6,0,t_res};
//right
Point(7) = {25+rad,0.6,0,t_res};
//top
Point(8) = {25,0.6+rad,0,t_res};
//bottom
Point(9) = {25,0.6-rad,0,t_res};

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