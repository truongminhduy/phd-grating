# DQNM expansion in diffraction grating

*Application of DQNM expansion for electromagnetic fields in the structure of diffraction grating.*
*Python and C++(for getdp and gmsh) code that reproduce the numerical results in the paper [doi:10.1109/COMPUMAG45669.2019.9032834](https://ieeexplore.ieee.org/document/9032834) and PhD manuscript [arXiv:2009.01307](https://arxiv.org/abs/2009.01307)*

## Affiliation

Minh Duy TRUONG: [minhduy.truong@fresnel.com](minhduy.truong@fresnel.com)
Aix Marseille Univ, CNRS, Centrale Marseille, Institut Fresnel, F-13013 Marseille, France
Athena team

## Instruction

1. Install getdp and gmsh and update the path leading to getdp and gmsh in the main files.
2. Run the main file "main_k00.py" and "main_k05.py"  in order to obtain data.
3. Run 6 files "plot_direct_1.py", "plot_eigen_improve.py", "plot_eigen_improve_2.py", "plot_energy_compare.py", "plot_direct_2.py" and "plot_compare_absorption.py" to plot the figures from data.

### Remark
* The data to plot the figures are saved in the .npz files which are named as follows: output_#angle#_#space2pml#a#pmlsize#_#expect_neig#.npz. 
* By changing these parameters in the main file "main_k00.py" such as #angle#, #space2pml#, #pmlsize#, #expect_neig# , it is possible to recover the data saved.
* The folder "save_k00_2" contains the eigenfrequency computation with no #EigFilter[]# for #kx0 = 0#.
* The folder "save_k00_1" contains the data of diffraction efficiency for #kx0 = 0#.
* The folder "save_k05" contains the data of diffraction efficiency for #kx0 = 0.5*pi/d#.
