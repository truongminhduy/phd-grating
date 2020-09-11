//////////////////////////////////////////////////////////////////////////////
Include "parameters_gmsh_getdp.dat";
Function{
  kx0 = -kx0;
}
Include "core_spectral.pro";

/////////// Post processing HZ EXY///////////////////////////////////////////////////////////////////
PostProcessing {
  { Name postpro_Hz; NameOfFormulation Hz;
    Quantity {
      { Name v ; Value { Local { [ {u}  ] ; In Omega; Jacobian JVol; } } }
      { Name dv; Value { Local { [ {Grad u}  ] ; In lamel_in; Jacobian JVol; } } }
      { Name EigenValuesReal; Value { Local{ [$EigenvalueReal]; In PrintPoint; Jacobian JVol; } } }
      { Name EigenValuesImag; Value { Local{ [$EigenvalueImag]; In PrintPoint; Jacobian JVol; } } }
    }
  }
  { Name postpro_Hz2; NameOfFormulation Hz2;
    Quantity {
      { Name EigenValuesReal; Value { Local{ [$EigenvalueReal]; In PrintPoint; Jacobian JVol; } } }
      { Name EigenValuesImag; Value { Local{ [$EigenvalueImag]; In PrintPoint; Jacobian JVol; } } }
    }
  }
}

PostOperation {
  { Name postop_Hz; NameOfPostProcessing postpro_Hz ;
    Operation {
      Print [EigenValuesReal, OnElementsOf PrintPoint, Format TimeTable, File "EigenValuesReal_ad.txt"];
      Print [EigenValuesImag, OnElementsOf PrintPoint, Format TimeTable, File "EigenValuesImag_ad.txt"];
      // Print [ v , OnElementsOf Omega, File "v.pos", EigenvalueLegend];
      For ki In {0:expect_neig-1}
        Print [ v , OnElementsOf Omega, TimeStep{ki}, File Sprintf["field/v_%g.pos",ki] ];
        Print [ dv, OnElementsOf lamel_in, TimeStep{ki}, File Sprintf["field/dv_%g.pos",ki] ];
      EndFor
    }
  }
  { Name postop_Hz2; NameOfPostProcessing postpro_Hz2 ;
    Operation {
      Print [EigenValuesReal, OnElementsOf PrintPoint, Format TimeTable, File "EigenValuesReal.txt"];
      Print [EigenValuesImag, OnElementsOf PrintPoint, Format TimeTable, File "EigenValuesImag.txt"];
    }
  }
}

