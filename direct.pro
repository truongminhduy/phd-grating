/////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
Include "core.pro";
////////////////////////////////////////////////////////////////////////

Formulation {
  {Name Hz_main; Type FemEquation;
    Quantity {
      { Name u; Type Local; NameOfSpace Hgrad;}
    }
    Equation {
      Galerkin { [ omega0^2*mur[]/cel^2                                                   * Dof{u}      , {u}     ] ; In Omega; Jacobian JVol; Integration Int_1;  }
      Galerkin { [-1/TensorDiag[ CompYY[epsilon[]] ,CompXX[epsilon[]] ,CompXX[epsilon[]] ]* Dof{Grad u} , {Grad u}] ; In Omega; Jacobian JVol; Integration Int_1; }                       
      Galerkin { [                                                                             source[] , {Grad u}] ; In Omega; Jacobian JVol; Integration Int_1;  }
    }
  }
}

Resolution {
  { Name Scattering;
    System {
      { Name S; NameOfFormulation Hz_main; Type ComplexValue; Frequency Freq;}
    }
    Operation {
      Generate[S];
      Solve[S];
      // SaveSolutions[S];
    }
  }
}

////// DATA PROCESS ////////////////////////////////////////////////////////////
PostProcessing {
  { Name get_Ed; NameOfFormulation Hz_main;
    Quantity {
      { Name us   ; Value { Local { [ {u}    ] ; In Omega; Jacobian JVol; } } }         
      { Name Hz_diff; Value { Local { [ {u}+u_1_d[] ]; In Omega; Jacobian JVol; } } }      
      { Name Q_lamel_in ; Value { Integral { [ 0.5 * epsilon0*omega0*Fabs[eps_lamel_in_im] * 
          (SquNorm[ CompY[{Grad u}]/(-I[]*omega0*epsilon0*CompXX[epsilon[]])+Ex_r[]/CompXX[epsilon[]]+Ex_i[]/CompXX[epsilon[]]] 
         + SquNorm[-CompX[{Grad u}]/(-I[]*omega0*epsilon0*CompXX[epsilon[]])+Ey_r[]/CompXX[epsilon[]]+Ey_i[]/CompXX[epsilon[]]] ) / 
          (Pinc*d) ] ; In lamel_in  ; Integration Int_1 ; Jacobian JVol ; } } }   
    }
  }
}

PostOperation {
  { Name postop_scat; NameOfPostProcessing get_Ed ;
    Operation {
      Print [ us, OnElementsOf Omega, File "us.pos" ];
      For i In {0:nb_slice-1}
        Print [Hz_diff , OnLine { {-d/2,ycut_sup~{i},0} {d/2, ycut_sup~{i},0} } {npt_integ-1}, File > "./Views/sup_field_cuts.out", Format Table];
        Print [Hz_diff , OnLine { {-d/2,ycut_sub~{i},0} {d/2, ycut_sub~{i},0} } {npt_integ-1}, File > "./Views/sub_field_cuts.out", Format Table];
      EndFor
      Print[ Q_lamel_in[lamel_in] , OnGlobal, File "./Views/temp-Q_lamel_in.txt"  , Format Table ];
    }
  }
}
///////////////////////////////////////////////////////////// END //////////////
////////////////////////////////////////////////////////////////////////////////
