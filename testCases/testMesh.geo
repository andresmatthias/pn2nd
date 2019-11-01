// run in command line:  gmsh testMesh.geo -2 
lc = 1;
Point(1) = {0, 0, 0, lc};
Point(2) = {1.0, 0, 0, lc};
Point(3) = {1.0, 1.0, 0, lc};
Point(4) = {0, 1.0, 0, lc};

Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 1};
Line(4) = {1, 2};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Physical Point(1) = {1, 2, 3, 4};
Physical Line("right",1) = {1}; 
Physical Line("top",2) = {2}; 
Physical Line("left",3) = {3};
Physical Line("bottom",4) = {4};
Physical Surface("2DLayer",999) = {1};
