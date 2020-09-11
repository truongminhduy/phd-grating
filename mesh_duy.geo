Include "parameters_gmsh_getdp.dat";

// fake_lambda = 400*nm;
fake_lambda = 20*nm;
// fake_lambda = 40*nm;
lc_lamel      = fake_lambda/Sqrt[epsilon_oo_re];
lc_lamel_out  = fake_lambda/Sqrt[eps_sup_re];
lc_layer_acc  = fake_lambda/Sqrt[eps_sub_re];
lc_layer_cov  = fake_lambda/Sqrt[eps_sup_re];
lc_sub		  = fake_lambda/Sqrt[eps_sub_re];
lc_sup		  = fake_lambda/Sqrt[eps_sup_re];
lc_pmlbot	  = fake_lambda/Sqrt[eps_sub_re];
lc_pmltop	  = fake_lambda/Sqrt[eps_sup_re];


Point(1)  = {-d/2., -space2pml-pmlsize/Sqrt[eps_sub_re]    , 0. , lc_pmlbot};
Point(2)  = {-d/2., -space2pml                             , 0. , lc_sub};
Point(3)  = {-d/2., 0.                                     , 0. , lc_sub};
Point(4)  = {-d/2., sub_acc                                , 0. , lc_layer_acc};
Point(5)  = {-d/2., sub_acc+d_y                            , 0. , lc_lamel_out};
Point(6)  = {-d/2., sub_acc+d_y+sup_cov                    , 0. , lc_layer_cov};
Point(7)  = {-d/2., sub_acc+d_y+sup_cov+space2pml          , 0. , lc_sup};
// Point(8)  = {-d/2., sub_acc+d_y+sup_cov+space2pml+pmlsize/Sqrt[eps_sup_re] , 0. , lc_pmltop};
Point(8)  = {-d/2., sub_acc+d_y+sup_cov+space2pml+pmlsize/Sqrt[eps_sub_re] , 0. , lc_pmltop};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};

out[]= Extrude{d,0,0}{Line{1};Line{2};Line{3};Line{4};Line{5};Line{6};Line{7};Layers{d/lc_sub};};

Delete{
	 Surface{11,15,19,23,27,31,35} ;
 	 Line{18,22};
      }

Point(40) = {-d_x/2, sub_acc     , 0 , lc_lamel};
Point(41) = { d_x/2, sub_acc     , 0 , lc_lamel};
Point(42) = {-d_x/2, sub_acc+ d_y, 0 , lc_lamel};
Point(43) = { d_x/2, sub_acc+ d_y, 0 , lc_lamel};

Line(50) = {4  , 40};
Line(51) = {40 , 41};
Line(52) = {41 , 14};
Line(53) = {5  , 42};
Line(54) = {42 , 43};
Line(55) = {43 , 16};
Line(56) = {40 , 42};
Line(57) = {41 , 43};

Line Loop(64) = {1, 10, -8, -9};
Plane Surface(65) = {-64};
Line Loop(66) = {2, 14, -12, -10};
Plane Surface(67) = {-66};
Line Loop(68) = {3, 50, 51, 52, -16, -14};
Plane Surface(69) = {-68};
Line Loop(70) = {4, 53, -56, -50};
Plane Surface(71) = {-70};
Line Loop(72) = {52, 20, -55, -57};
Plane Surface(73) = {72};
Line Loop(74) = {51, 57, -54, -56};
Plane Surface(75) = {74};
Line Loop(80) = {5, 26, -24, -55, -54, -53};
Plane Surface(81) = {-80};
Line Loop(82) = {6, 30, -28, -26};
Plane Surface(83) = {-82};
Line Loop(84) = {7, 34, -32, -30};
Plane Surface(85) = {-84};

Physical Line(101) = {1, 2, 3, 4, 5, 6, 7};		        // Bloch_LeftX-
Physical Line(102) = {8,12,16,20,24,28,32};		        // Bloch_RightX+
Physical Line(110) = {9, 34};			                // Dirichlet

Physical Surface(1000)  = {65};        // pmlbot
Physical Surface(2000)  = {67};        // sub
Physical Surface(3000)  = {69};        // layer_acc
Physical Surface(4000)  = {71,73};     // lamel_out_left
Physical Surface(5000)  = {75};  // lamel_in
Physical Surface(6000)  = {81};        // layer_cov
Physical Surface(7000) = {83};        // sup
Physical Surface(8000) = {85};        // pmltop
Physical Point(10000) = {7};

Point(51)  = {-d/2. ,ycut_sup_max, 0. , lc_sub};
Point(52)  = {-d/2. ,ycut_sup_min, 0. , lc_sub};
Point(53)  = {-d/2. ,ycut_sub_max, 0. , lc_sub};
Point(54)  = {-d/2. ,ycut_sub_min, 0. , lc_sub};
Physical Point(10001) = {51};  
Physical Point(10002) = {52};  
Physical Point(10003) = {53};  
Physical Point(10004) = {54}; 