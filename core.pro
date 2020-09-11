/////////////////////////////////////////////////////////////////////////
Group {
  pmlbot   = Region[1001];
  sub      = Region[1002];
  lamel_in = Region[1003];
  sup      = Region[1004];
  pmltop   = Region[1005];
  Omega_out = Region[{pmltop,pmlbot,sup,sub}];
  Omega     = Region[{Omega_out,lamel_in}];

  SurfBlochLeft  = Region[101];
  SurfBlochRight = Region[102];
  SurfDirichlet  = Region[110];
  PrintPoint     = Region[2000];
}
/////////////////////////////////////////////////////////////////////////////
Function{
    I[]    = Complex[0.0,1.0];
    Freq   = cel/lambda0;
    omega0 = 2.*Pi*cel/lambda0;
    k0     = 2.*Pi/lambda0;        
    k_sup  = k0*Sqrt[eps_sup_re];
    k_sub  = k0*Sqrt[eps_sub_re];

    alpha_sup = kx0;
    theta     = Asin[kx0/k_sup];
    beta_sup  = k_sup*Cos[theta];
    beta_sub  = Sqrt[k_sub*k_sub-alpha_sup*alpha_sup];

    beta_S_sup[] = beta_sup/Complex[eps_sup_re,eps_sup_im];
    beta_S_sub[] = beta_sub/Complex[eps_sub_re,eps_sub_im];
    R[]          = (beta_S_sup[]-beta_S_sub[])/(beta_S_sup[]+beta_S_sub[]);
    T[]          = (2.*beta_S_sup[])/(beta_S_sup[]+beta_S_sub[]);
    Pinc         = 0.5*A*A*Sqrt[mu0/epsilon0] * Cos[theta]; // total energy

    deph[] = Complex[ Cos[alpha_sup*d] , Sin[alpha_sup*d] ];
    // sy[]   = Complex[Cos[angle],Sin[angle]];
    sy[]   = Complex[1,Sin[angle]];

    epsilon[pmltop]   = TensorDiag[ eps_sup_re*sy[], eps_sup_re/sy[], eps_sup_re*sy[] ];
    epsilon[sup]      = Complex[eps_sup_re,eps_sup_im] * TensorDiag[1,1,1];
    epsilon[lamel_in] = Complex[eps_lamel_in_re,eps_lamel_in_im] * TensorDiag[1,1,1];
    epsilon[sub]      = Complex[eps_sub_re,eps_sub_im] * TensorDiag[1,1,1];
    epsilon[pmlbot]   = TensorDiag[ eps_sub_re*sy[], eps_sub_re/sy[], eps_sub_re*sy[] ];

    epsilon_annex[pmltop]   = TensorDiag[ eps_sup_re*sy[], eps_sup_re/sy[], eps_sup_re*sy[] ];
    epsilon_annex[sup]      = Complex[eps_sup_re,eps_sup_im] * TensorDiag[1,1,1];
    epsilon_annex[lamel_in] = Complex[eps_sup_re,eps_sup_im] * TensorDiag[1,1,1];
    epsilon_annex[sub]      = Complex[eps_sub_re,eps_sub_im] * TensorDiag[1,1,1];
    epsilon_annex[pmlbot]   = TensorDiag[ eps_sub_re*sy[], eps_sub_re/sy[], eps_sub_re*sy[] ];
    
    epsilon1[pmltop]   = TensorDiag[ eps_sup_re*sy[], eps_sup_re/sy[], eps_sup_re*sy[] ];
    epsilon1[sup]      = TensorDiag[1,1,1];
    epsilon1[lamel_in] = TensorDiag[1,1,1];  
    epsilon1[sub]      = TensorDiag[1,1,1];
    epsilon1[pmlbot]   = TensorDiag[ eps_sub_re*sy[], eps_sub_re/sy[], eps_sub_re*sy[] ];

    mur[pmltop]   = TensorDiag[sy[],1/sy[],sy[]];
    mur[sup]      = TensorDiag[1,1,1];
    mur[lamel_in] = TensorDiag[1,1,1];
    mur[sub]      = TensorDiag[1,1,1];
    mur[pmlbot]   = TensorDiag[sy[],1/sy[],sy[]];

    // u is E electric field in this case ///////////////////////////////////////
    // incident field
    u_i[pmltop]   = 0.;
    u_i[sup]      = A*Complex[ Cos[alpha_sup*X[]-beta_sup*Y[]] , Sin[alpha_sup*X[]-beta_sup*Y[]] ];
    u_i[lamel_in] = A*Complex[ Cos[alpha_sup*X[]-beta_sup*Y[]] , Sin[alpha_sup*X[]-beta_sup*Y[]] ];
    u_i[sub]      = 0.;
    u_i[pmlbot]   = 0.;

    // refected field
    u_r[pmltop]   = 0.;
    u_r[sup]      = R[]*Complex[ Cos[alpha_sup*X[]+beta_sup*Y[]] , Sin[alpha_sup*X[]+beta_sup*Y[]] ];
    u_r[lamel_in] = R[]*Complex[ Cos[alpha_sup*X[]+beta_sup*Y[]] , Sin[alpha_sup*X[]+beta_sup*Y[]] ];
    u_r[sub]      = 0;
    u_r[pmlbot]   = 0.;

    // tranmitted field
    u_t[pmltop]   = 0.;
    u_t[sup]      = 0.;
    u_t[lamel_in] = 0.;
    u_t[sub]      = T[]*Complex[ Cos[alpha_sup*X[]-beta_sub*Y[]] , Sin[alpha_sup*X[]-beta_sub*Y[]] ];
    u_t[pmlbot]   = 0.;

    u_1[]    = u_i[]+u_r[]+u_t[];
    u_1_d[]  = u_r[]+u_t[];

    //  from u = H we compute E
    Ex_i[]  =  (-I[]*beta_sup *u_i[]) / (-I[]*omega0*epsilon0*CompXX[epsilon_annex[]]);
    Ey_i[]  = -( I[]*alpha_sup*u_i[]) / (-I[]*omega0*epsilon0*CompXX[epsilon_annex[]]);
    Ex_r[]  =  ( I[]*beta_sup *u_r[]) / (-I[]*omega0*epsilon0*CompXX[epsilon_annex[]]);
    Ey_r[]  = -( I[]*alpha_sup*u_r[]) / (-I[]*omega0*epsilon0*CompXX[epsilon_annex[]]);
    Ex_t[]  =  (-I[]*beta_sub *u_t[]) / (-I[]*omega0*epsilon0*CompXX[epsilon_annex[]]);
    Ey_t[]  = -( I[]*alpha_sup*u_t[]) / (-I[]*omega0*epsilon0*CompXX[epsilon_annex[]]);

    source_i[] = TensorDiag[alpha_sup,-beta_sup,0.]* u_i[];       
    source_r[] = TensorDiag[alpha_sup, beta_sup,0.]* u_r[];
    source[]   = I[] *(1/epsilon_annex[]-1/epsilon[]) *(source_i[]+source_r[]);
    
    E_i[] = TensorDiag[Ex_i[],Ey_i[],0];
    E_r[] = TensorDiag[Ex_r[],Ey_r[],0];
    E_t[] = TensorDiag[Ex_t[],Ey_t[],0];
    EE_i[] = Vector[Ex_i[],Ey_i[],0];
    EE[]  = Vector[ Ex_i[]+Ex_r[]+Ex_t[] , Ey_i[]+Ey_r[]+Ey_t[] , 0 ];
    sourceE[] = (omega0/cel)^2*(epsilon_annex[]-epsilon[])*(E_i[]+E_r[]);
    
    // SE[] = Vector[CompXX[sourceE[]],CompYY[sourceE[]],0];
    SE[] = Vector[(omega0/cel)^2*(CompZZ[epsilon_annex[]]-CompZZ[epsilon[]])*(Ex_i[]+Ex_r[])
                 ,(omega0/cel)^2*(CompZZ[epsilon_annex[]]-CompZZ[epsilon[]])*(Ey_i[]+Ey_r[]),0];

    // y-position to compute the diffraction coefficients
    For i In {0:nb_slice-1}     // loop over all the slides (we take the mean value of all diff coefficients wrt to y-position later)
        ycut_sub~{i} = ycut_sub_min + i*(ycut_sub_max-ycut_sub_min)/(nb_slice-1);
        ycut_sup~{i} = ycut_sup_min + i*(ycut_sup_max-ycut_sup_min)/(nb_slice-1);
    EndFor
}

//// FUNCTION //////////////////////////////////////////////////////////////////
Constraint {
    {Name Dirichlet; Type Assign;
        Case {
            { Region SurfDirichlet; Value 0.; }
        }
    }
    {Name BlochX;
        Case {
            { Region SurfBlochRight; Type LinkCplx ; RegionRef SurfBlochLeft; Coefficient deph[]; Function Vector[$X-d,$Y,$Z] ; }
        }
    }
}

Jacobian {
  { Name JVol ; Case { { Region All ; Jacobian Vol ; } } }
  { Name JSur ; Case { { Region All ; Jacobian Sur ; } } }
  { Name JLin ; Case { { Region All ; Jacobian Lin ; } } }
}

Integration {
    { Name Int_1 ;
        Case {
            { Type Gauss ;
                Case {
                    { GeoElement Point   ; NumberOfPoints  1 ; }
                    { GeoElement Line    ; NumberOfPoints  4 ; }
                    { GeoElement Triangle; NumberOfPoints  12 ; }
                }
            }
        }
    }
}

FunctionSpace {
  { Name Hgrad; Type Form0;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Node   ; Support Region[Omega]; Entity NodesOf[All]; }
      { Name sn2; NameOfCoef un2; Function BF_Node_2E; Support Region[Omega]; Entity EdgesOf[All]; }
    }
    Constraint {
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint BlochX; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint BlochX; }
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
    }
  }
    { Name Hgrad_perp; Type Form1P;
    BasisFunction {
      { Name un;  NameOfCoef un;  Function BF_PerpendicularEdge_1N; Support Region[Omega]; Entity NodesOf[Omega]; }
      { Name un2; NameOfCoef un2; Function BF_PerpendicularEdge_2E; Support Region[Omega]; Entity EdgesOf[Omega]; }
     }
    Constraint {
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint BlochX; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint BlochX; }
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
    }
  }

  { Name Hvector; Type Form1P;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_PerpendicularEdge_1N; Support Region[Omega]; Entity NodesOf[All]; }
      { Name sn2; NameOfCoef un2; Function BF_PerpendicularEdge_2E; Support Region[Omega]; Entity EdgesOf[All]; }
    }
    Constraint {
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint BlochX; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint BlochX; }
      { NameOfCoef un;  EntityType NodesOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef un2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
    }
  }

 { Name E_Edge; Type Form1;
    BasisFunction {
      { Name sE;   NameOfCoef uE;  Function BF_Edge;    Support Region[Omega]; Entity EdgesOf[Omega]; }
      { Name sE2;  NameOfCoef uE2; Function BF_Edge_2E; Support Region[Omega]; Entity EdgesOf[Omega]; }
  }
  Constraint {
      { NameOfCoef uE;  EntityType EdgesOf ; NameOfConstraint BlochX; }
      { NameOfCoef uE2; EntityType EdgesOf ; NameOfConstraint BlochX; }
      { NameOfCoef uE;  EntityType EdgesOf ; NameOfConstraint Dirichlet; }
      { NameOfCoef uE2; EntityType EdgesOf ; NameOfConstraint Dirichlet; }
    }
  }

  { Name EaD1_Edge; Type Form1;
    BasisFunction {
      { Name sEaD1;   NameOfCoef uEaD1;  Function BF_Edge;    Support Region[lamel_in]; Entity EdgesOf[lamel_in]; }
      { Name sEaD12;  NameOfCoef uEaD12; Function BF_Edge_2E; Support Region[lamel_in]; Entity EdgesOf[lamel_in]; }
    }
  }
}
