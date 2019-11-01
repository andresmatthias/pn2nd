lc = 0.1;
Point(1) = {0, 0, 0, lc};
Point(2) = {3.0, 0, 0, lc};
Point(3) = {3.0, 1.0, 0, lc};
Point(4) = {0, 1.0, 0, lc};

Point(5) = {1.0, 0.0, 0, lc};
Point(6) = {1.0, 1.0, 0, lc};

Point(7) = {0.45, 0, 0, lc};
Point(8) = {0.85, 0, 0, lc};
Point(9) = {0.45, 0.55, 0, lc};
Point(10) = {0.85, 0.55, 0, lc};

Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 1};
Line(4) = {1, 2};

Line(5) = {1, 5};
Line(6) = {4, 6};
Line(7) = {5, 6};

Line(8) = {7, 9};
Line(9) = {9, 10};
Line(10) = {8, 10};


Line Loop(1) = {1, 2, 3, 4};

// attractors
Field[1] = Attractor;
Field[1].NNodesByEdge = 100;
Field[1].EdgesList = {3, 5, 6};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc / 2;
Field[2].LcMax = lc;
Field[2].DistMin = 0.0;
Field[2].DistMax = 2;

Field[3] = Attractor;
Field[3].NNodesByEdge = 100;
Field[3].EdgesList = {8, 9, 10};
Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = lc / 10;
Field[4].LcMax = lc;
Field[4].DistMin = 0.0;
Field[4].DistMax = 0.5;

Field[5] = Min;
Field[5].FieldsList = {2, 4};

Background Field = 5;

Plane Surface(1) = {1};

Physical Point(1) = {1, 2, 3, 4};
Physical Line("right",1) = {1}; 
Physical Line("top",2) = {2}; 
Physical Line("left",3) = {3};
Physical Line("bottom",4) = {4};
Physical Surface("2DLayer",999) = {1};

