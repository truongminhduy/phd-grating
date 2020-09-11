
Plugin(MathEval).Expression0 = "Sqrt((v0-w0)^2)"; //+(v1-w1)^2
Plugin(MathEval).TimeStep = 0;
Plugin(MathEval).View = 0;
Plugin(MathEval).OtherTimeStep = 0;
Plugin(MathEval).OtherView = 1;
Plugin(MathEval).Run;
ve = PostProcessing.NbViews-1;
View[ve].Name = "error on modal estimation";
