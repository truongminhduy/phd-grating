//////////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
Include "core_spectral.pro";
///////////////////////////////////////////////////////////////////////////////////////////////

PostProcessing {
  { Name postpro_norm; NameOfFormulation Hz;
    Quantity {
      { Name norm1 ; Value { Integral { [ CompZZ[mur[]]* {u}*{u} ]; In Omega; Integration Int_1; Jacobian JVol; } } }   
      { Name norm2 ; Value { Integral { [ {Grad u}*{Grad u} ] ; In lamel_in; Integration Int_1; Jacobian JVol; } } } 
    }
  }
}

PostOperation {
  { Name postop_norm; NameOfPostProcessing postpro_norm ;
    Operation {
      Print [  norm1[Omega]    , OnElementsOf PrintPoint, Format TimeTable, File "norm1.txt"];
      Print [  norm2[lamel_in] , OnElementsOf PrintPoint, Format TimeTable, File "norm2.txt"];
    }
  }
}
