//////////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
Include "Pns.dat";
Include "core.pro";
Function{
  For n In {0:neig-1}
    P~{n}[] = Complex[ Pns_re~{n}, Pns_im~{n} ];
  EndFor
}

////////////////////////////////////////////////////////////////////////////////////////////////
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
}

Resolution {
  { Name Projection;
    System{
      { Name M1; NameOfFormulation Hz ; Type ComplexValue; }
    }
    Operation{
      GenerateSeparate[M1];
      EigenSolve[M1,neig,eig_target_re,eig_target_im];  
      SaveSolutions[M1];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
PostProcessing {
  { Name postpro_modal; NameOfFormulation Hz;
    Quantity {
      { Name up;
        Value {
          For i In {0:neig-1}
            Local { [ P~{i}[] * {u}[neig-1-i] ]; In Omega; Jacobian JVol; }
            // Local { [ P~{i}[] * CompZ[{u}[neig-1-i]] ]; In Omega; Jacobian JVol; }
          EndFor

          // Local { [ P~{22}[] * {u}[neig-1-22] ]; In Omega; Jacobian JVol; }
        }
      }
    }
  }
}

PostOperation {
  { Name postop_modal; NameOfPostProcessing postpro_modal ;
    Operation {
      Print [ up, OnElementsOf Omega, File "up.pos", LastTimeStepOnly ];
    }
  }
}
