//////////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////
Function{
  I[] = Complex[0,1];
  deph[] = Complex[ Cos[kx0*d] , Sin[kx0*d] ];

  //////// SPECTRAL PROBLEM ///////////////////////////////////////////////////////////////////
  sy[]  = Complex[Cos[angle],Sin[angle]];

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

  EigFilter[] = (Norm[$EigenvalueReal] > 1e-4);
}

Constraint {
  {Name BlochX;
    Case { { Region SurfBlochRight; Type LinkCplx ; RegionRef SurfBlochLeft; Coefficient deph[]; Function Vector[$X-d,$Y,$Z] ; }}
  }
  {Name Dirichlet; Type Assign;
    Case { { Region SurfDirichlet ; Value 0.; }}
  }
}

Jacobian {
  { Name JVol ; Case { { Region All ; Jacobian Vol ; } } }
}

Integration {
  { Name Int_1 ;
    Case {
      { Type Gauss ;
        Case {
          { GeoElement Point       ; NumberOfPoints  1 ; }
          { GeoElement Line        ; NumberOfPoints  4 ; }
          { GeoElement Triangle    ; NumberOfPoints  6 ; }
        }
      }
    }
  }
}

/////// DEFINE FUNCTIONSPACE AND WEAK EQUATION /////////////////////////////////
FunctionSpace {
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
  { Name Hgrad; Type Form0;
    BasisFunction {
      { Name sn;  NameOfCoef un;  Function BF_Node   ; Support Region[Omega]; Entity NodesOf[Omega]; }
      { Name sn2; NameOfCoef un2; Function BF_Node_2E; Support Region[Omega]; Entity EdgesOf[Omega]; }
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
}

////////////////////////////////////////////////////////////////////////////////////////////
Formulation {
  // {Name Hz; Type FemEquation; // just the name dontcare about it
  //   Quantity {{ Name u   ; Type Local; NameOfSpace Hgrad_perp  ;}}
  //   Equation {
  //     Galerkin {    [ 1/epsilon1[] * omega_p^2             * Dof{Curl u}, {Curl u} ];          In Omega_out; Jacobian JVol; Integration Int_1;}
  //     Galerkin { Eig[-1/epsilon1[] * gamma*  epsilon_oo_re * Dof{Curl u}, {Curl u} ]; Order 1; In Omega_out; Jacobian JVol; Integration Int_1;}
  //     Galerkin { Eig[ 1/epsilon1[] *         epsilon_oo_re * Dof{Curl u}, {Curl u} ]; Order 2; In Omega_out; Jacobian JVol; Integration Int_1;}
  //     Galerkin { Eig[-1/epsilon1[] * gamma                 * Dof{Curl u}, {Curl u} ]; Order 1; In lamel_in; Jacobian JVol; Integration Int_1;}
  //     Galerkin { Eig[ 1/epsilon1[]                         * Dof{Curl u}, {Curl u} ]; Order 2; In lamel_in; Jacobian JVol; Integration Int_1;}
  //     Galerkin { Eig[ mur[]/cel^2   * omega_p^2             *      Dof{u},      {u} ]; Order 2; In Omega; Jacobian JVol; Integration Int_1;}
  //     Galerkin { Eig[-mur[]/cel^2   * gamma*  epsilon_oo_re *      Dof{u},      {u} ]; Order 3; In Omega; Jacobian JVol; Integration Int_1;}
  //     Galerkin { Eig[ mur[]/cel^2   *         epsilon_oo_re *      Dof{u},      {u} ]; Order 4; In Omega; Jacobian JVol; Integration Int_1;}
  //   }
  // }
  {Name Hz; Type FemEquation;
    Quantity {{ Name u   ; Type Local; NameOfSpace Hgrad  ;}}
    Equation {
      Galerkin {    [ 1/TensorDiag[ CompYY[epsilon1[]] ,CompXX[epsilon1[]] ,CompXX[epsilon1[]] ] * omega_p^2             * Dof{Grad u}, {Grad u} ];          In Omega_out; Jacobian JVol; Integration Int_1;}
      Galerkin { Eig[-1/TensorDiag[ CompYY[epsilon1[]] ,CompXX[epsilon1[]] ,CompXX[epsilon1[]] ] * gamma*  epsilon_oo_re * Dof{Grad u}, {Grad u} ]; Order 1; In Omega_out; Jacobian JVol; Integration Int_1;}
      Galerkin { Eig[ 1/TensorDiag[ CompYY[epsilon1[]] ,CompXX[epsilon1[]] ,CompXX[epsilon1[]] ] *         epsilon_oo_re * Dof{Grad u}, {Grad u} ]; Order 2; In Omega_out; Jacobian JVol; Integration Int_1;}
      Galerkin { Eig[-1/TensorDiag[ CompYY[epsilon1[]] ,CompXX[epsilon1[]] ,CompXX[epsilon1[]] ] * gamma                 * Dof{Grad u}, {Grad u} ]; Order 1; In lamel_in; Jacobian JVol; Integration Int_1;}
      Galerkin { Eig[ 1/TensorDiag[ CompYY[epsilon1[]] ,CompXX[epsilon1[]] ,CompXX[epsilon1[]] ]                         * Dof{Grad u}, {Grad u} ]; Order 2; In lamel_in; Jacobian JVol; Integration Int_1;}
      Galerkin { Eig[ mur[]/cel^2   * omega_p^2             *      Dof{u},      {u} ]; Order 2; In Omega; Jacobian JVol; Integration Int_1;}
      Galerkin { Eig[-mur[]/cel^2   * gamma*  epsilon_oo_re *      Dof{u},      {u} ]; Order 3; In Omega; Jacobian JVol; Integration Int_1;}
      Galerkin { Eig[ mur[]/cel^2   *         epsilon_oo_re *      Dof{u},      {u} ]; Order 4; In Omega; Jacobian JVol; Integration Int_1;}
    }
  }
  {Name Hz2; Type FemEquation; 
    Quantity {{ Name u   ; Type Local; NameOfSpace Hgrad_perp  ;}}
    Equation {
      Galerkin { Eig[ 1/epsilon1[] * Dof{Curl u}, {Curl u} ]; Rational 1; In Omega_out; Jacobian JVol; Integration Int_1;}
      Galerkin { Eig[ 1/epsilon1[] * Dof{Curl u}, {Curl u} ]; Rational 2; In lamel_in ; Jacobian JVol; Integration Int_1;}
      Galerkin { Eig[ mur[]/cel^2  * Dof{u},      {u}      ]; Rational 3; In Omega    ; Jacobian JVol; Integration Int_1;}
    }
  }
}

Resolution {
  { Name Projection;
    System{
      { Name M1; NameOfFormulation Hz ; Type ComplexValue; }
    }
    Operation{
      For i In {0:neig-1}
        GmshRead[Sprintf["field/u_%g.pos",i]];
        GmshRead[Sprintf["field/v_%g.pos",i]];
        GmshRead[Sprintf["field/du_%g.pos",i]];
        GmshRead[Sprintf["field/dv_%g.pos",i]];
      EndFor
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
PostProcessing {
  { Name postpro_norm; NameOfFormulation Hz;
    Quantity {
      For i In {0:neig-1}
        { Name norma~{i} ; Value { Integral { [ CompZZ[mur[]]*Field[XYZ[]]{4*i+0}*Field[XYZ[]]{4*i+1} ]; In Omega; Integration Int_1; Jacobian JVol; } } } 
        { Name normb~{i} ; Value { Integral { [ Field[XYZ[]]{4*i+2}*Field[XYZ[]]{4*i+3} ]; In lamel_in; Integration Int_1; Jacobian JVol; } } }    
      EndFor
    }
  }
}

PostOperation {
  { Name postop_norm; NameOfPostProcessing postpro_norm ;
    Operation {
      For i In {0:neig-1}
        Print [  norma~{i}[Omega]   , OnElementsOf PrintPoint, Format TimeTable, File > "norm1.txt"];
        Print [  normb~{i}[lamel_in], OnElementsOf PrintPoint, Format TimeTable, File > "norm2.txt"];
      EndFor
    }
  }
}
