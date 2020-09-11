Include "parameters_gmsh_getdp.dat";

// fake_lambda   = 6e-2;
fake_lambda   = 3e-2;
lc_lamel      = 0.8*fake_lambda/Sqrt[epsilon_oo_re];
lc_layer_acc  =     fake_lambda/Sqrt[eps_sub_re];
lc_sub		  =     fake_lambda/Sqrt[eps_sub_re];
lc_sup		  =     fake_lambda/Sqrt[eps_sup_re];
// lc_pmlbot	  =   4*fake_lambda/Sqrt[eps_sub_re];
// lc_pmltop	  =   4*fake_lambda/Sqrt[eps_sup_re];
lc_pmlbot	  = 1.5*fake_lambda/Sqrt[eps_sub_re];
lc_pmltop	  = 1.5*fake_lambda/Sqrt[eps_sup_re];


Point(1)  = {-d/2.,     -space2pml -pmlsize/Sqrt[eps_sub_re], 0. , lc_pmlbot};
Point(2)  = {-d/2.,     -space2pml                          , 0. , lc_sub};
Point(3)  = {-d/2., 0                                       , 0. , lc_layer_acc};
Point(4)  = {-d/2., d_y +space2pml                          , 0. , lc_sup};
Point(5)  = {-d/2., d_y +space2pml +pmlsize/Sqrt[eps_sub_re], 0. , lc_pmltop};

Point(6)  = { d/2.,     -space2pml -pmlsize/Sqrt[eps_sub_re], 0. , lc_pmlbot};
Point(7)  = { d/2.,     -space2pml                          , 0. , lc_sub};
Point(8)  = { d/2., 0                                       , 0. , lc_layer_acc};
Point(9)  = { d/2., d_y +space2pml                          , 0. , lc_sup};
Point(10) = { d/2., d_y +space2pml +pmlsize/Sqrt[eps_sub_re], 0. , lc_pmltop};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};

Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9,10};

Line(11) = {1, 6};
Line(12) = {2, 7};
Line(13) = {4, 9};
Line(14) = {5,10};

Point(11) = {-d_x/2, 0  , 0 , lc_lamel};
Point(12) = { d_x/2, 0  , 0 , lc_lamel};
Point(13) = {-d_x/2, d_y, 0 , lc_lamel};
Point(14) = { d_x/2, d_y, 0 , lc_lamel};

Line(20) = {3  , 11};
Line(21) = {11 , 12};
Line(22) = {12 ,  8};
Line(23) = {11 , 13};
Line(24) = {13 , 14};
Line(25) = {14 , 12};

Line Loop(30) = {-1, 11, 5, -12};
Line Loop(32) = {-2, 12, 6, -20, -21 , -22};
Line Loop(34) = {21, -23, -24, -25};
Line Loop(36) = {-3, 20, 22, 23, 24, 25, 7, -13};
Line Loop(38) = {-4, 13, 8, -14};

Plane Surface(31) = {30};
Plane Surface(33) = {32};
Plane Surface(35) = {34};
Plane Surface(37) = {36};
Plane Surface(39) = {38};


Physical Line(101) = {1, 2, 3, 4};	// Bloch_LeftX-
Physical Line(102) = {5, 6, 7, 8};	// Bloch_RightX+
Physical Line(110) = {11, 14};		// Dirichlet

Physical Surface(1001) = {31};  // pmlbot
Physical Surface(1002) = {33};  // sub
Physical Surface(1003) = {35};  // lamel_in
Physical Surface(1004) = {37};  // sup
Physical Surface(1005) = {39};  // pmltop
Physical Point(2000) = {1};