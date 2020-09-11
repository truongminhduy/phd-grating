+ Install getdp and gmsh and update the path leading to getdp and gmsh in the main*.py files.
+ Run the main file "main_k00.py" and "main_k05.py"  in order to obtain data.
+ Run 6 files "plot_direct_1.py", "plot_eigen_improve.py", "plot_eigen_improve_2.py", "plot_energy_compare.py", "plot_direct_2.py" and "plot_compare_absorption.py" to plot the figures from data.

- The data to plot the figures are saved in the .npz files which are named as follows: output_#angle#_#space2pml#a#pmlsize#_#expect_neig#.npz. 
- By changing these parameters in the main file "main_k00.py" such as #angle#, #space2pml#, #pmlsize#, #expect_neig# , it is possible to recover the data saved.
- The folder "save_k00_2" contains the eigenfrequency computation with no #EigFilter[]# for #kx0 = 0#.
- The folder "save_k00_1" contains the data of diffraction efficiency for #kx0 = 0#.
- The folder "save_k05" contains the data of diffraction efficiency for #kx0 = 0.5*pi/d#.
