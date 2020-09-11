//////////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
Include "core.pro";
////////////////////////////////////////////////////////////////////////////////////////////////
Formulation {
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
  // {Name Hz2; Type FemEquation; 
  //   Quantity {{ Name u   ; Type Local; NameOfSpace Hgrad_perp  ;}}
  //   Equation {
  //     Galerkin { Eig[ 1/epsilon1[] * Dof{Curl u}, {Curl u} ]; Rational 1; In Omega_out; Jacobian JVol; Integration Int_1;}
  //     Galerkin { Eig[ 1/epsilon1[] * Dof{Curl u}, {Curl u} ]; Rational 2; In lamel_in ; Jacobian JVol; Integration Int_1;}
  //     Galerkin { Eig[ mur[]/cel^2  * Dof{u},      {u}      ]; Rational 3; In Omega    ; Jacobian JVol; Integration Int_1;}
  //   }
  // }
}

Resolution {
  { Name Projection;
    System{
      { Name M1; NameOfFormulation Hz; Type ComplexValue; }
      // { Name M1; NameOfFormulation Hz; Type ComplexValue; }
    }
    Operation{
      GenerateSeparate[M1];
      // EigenSolve[M1,expect_neig,eig_target_re,eig_target_im,EigFilter[]];  
      EigenSolve[M1,expect_neig,eig_target_re,eig_target_im];  
      SaveSolutions[M1];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
PostProcessing {
  { Name postpro_Jn; NameOfFormulation Hz;
    Quantity {
      { Name Jn ; Value { Integral { [ - ( CompXX[source[]]*CompX[{Grad u}] + CompYY[source[]]*CompY[{Grad u}] ) ]; In Omega; Integration Int_1; Jacobian JVol; } } }
      // { Name Jn ; Value { Integral { [ CompXX[source[]]*CompY[{Curl u}] - CompYY[source[]]*CompX[{Curl u}]  ]; In Omega; Integration Int_1; Jacobian JVol; } } }
    }
  }
}

PostOperation {
  { Name postop_Jn; NameOfPostProcessing postpro_Jn ;
    Operation {
      Print [  Jn[Omega], OnElementsOf PrintPoint, Format TimeTable, File "Jns.txt"];
    }
  }
}
