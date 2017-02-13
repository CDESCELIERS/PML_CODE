lc = 0.3/(80.0/24.0);
lc2 = 0.5/(80.0/24.0);

Lright    = 12 ;
Lbottom   = 24;
y_PML     = -10;
y_layer_1 = 0;
y_layer_2 = -0.15;
y_layer_3 = -2;
y_layer_4 = -5;

Point(1)  = {-Lright, y_layer_1, 0, lc};
Point(2)  = { Lright, y_layer_1, 0, lc};
Point(3)  = { Lright, y_layer_2, 0, lc};
Point(4)  = {-Lright, y_layer_2, 0, lc};
Point(5)  = {-Lright, y_layer_3, 0, lc};
Point(6)  = { Lright, y_layer_3, 0, lc};
Point(7)  = {-Lright, y_layer_4, 0, lc2};
Point(8)  = { Lright, y_layer_4, 0, lc2};
Point(9)  = {-Lright, y_PML    , 0, lc2};
Point(10) = { Lright, y_PML    , 0, lc2};
Point(11) = {-Lright, -Lbottom , 0, lc2};
Point(12) = { Lright, -Lbottom , 0, lc2};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line(7) = {6, 3};
Line(8) = {6, 8};
Line(9) = {8, 7};
Line(10) = {7, 5};
Line(11) = {7, 9};
Line(12) = {9, 10};
Line(13) = {10, 8};
Line(14) = {10, 12};
Line(15) = {12, 11};
Line(16) = {11, 9};

Line Loop(17) = {4, 1, 2, 3};
Ruled Surface(18) = {17};

Line Loop(19) = {7, 3, 5, 6};
Ruled Surface(20) = {19};

Line Loop(21) = {10, 6, 8, 9};
Ruled Surface(22) = {21};

Line Loop(23) = {13, 9, 11, 12};
Ruled Surface(24) = {23};

Line Loop(25) = {16, 12, 14, 15};
Ruled Surface(26) = {25};



