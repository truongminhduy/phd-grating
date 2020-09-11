################################################################################
#### GRATING PROBLEM ###########################################################
####################################################################!!!#########
import numpy as np
import scipy as sc
import scipy.io as scio
import pylab as pl
import os
import sys

# from direct_diffraction_efficiencies_TE import *
from direct_diffraction_efficiencies_TM import *
import matplotlib.pyplot as pl
pi = np.pi

#### DEFINE FUNCTION ###########################################################
def dumpglob2npz(filename, current_globals):
    dic2save = {}
    for glob_key in sorted(current_globals.keys()):  
        glob_type = type(current_globals[glob_key])
        if glob_type in (int, bool, float, np.float64, complex, np.ndarray, str, list, dict) and glob_key[0]!='_':
            dic2save[glob_key] = current_globals[glob_key]
    np.savez_compressed(filename,saveddic = dic2save)
##################################################################################

#### MAIN 1 #######################################################################
### define inpute files ###########################################################
# str_gmsh_path  = '../../../Downloads/gmsh-4.3.0-Linux64/bin/'
# str_getdp_path = '../../../Downloads/getdp-3.1.0-Linux64/bin/'
str_gmsh_path  = ''
str_getdp_path = ''
mesh_filename = 'mesh_duy_simple.msh'
mesh_geo      = 'mesh_duy_simple.geo'
#### input #########################################################################
cc    = 299792458
mu_0  = pi*4e-7
eps_0 = 1./(mu_0*cc**2)
A     = 1.
scale_light = 1e8
micro       = 1e-6
scale_frequ = scale_light/micro # =1e14
cel      = cc/scale_light
mu0      = mu_0  *scale_light
epsilon0 = eps_0 *scale_light

N_d_order = 2
d         = 0.4825
d_y       = 0.13
d_x       = 0.3475          
space2pml = 2
pmlsize   = 25
# space2pml = 8.8
# pmlsize   = 11
### silver ####
epsilon_oo = 3.36174
omega_p    = 1.3388e16/scale_frequ
gamma      = 7.07592e13/scale_frequ
eps_sup    = 1
eps_sub    = 1  

scan_dist    = 1e-2
ycut_sup_max = d_y+space2pml-scan_dist
ycut_sup_min = d_y          +scan_dist
ycut_sub_max = -scan_dist
ycut_sub_min = -space2pml+scan_dist
npt_integ    = 100
nb_slice     = 5
expect_neig = 720
# expect_neig = 300
eig_target_im = 2.*pi*cel/(0.6)
# eig_target_im = 2.*pi*cel/(0.7)
eig_target_re = 0
###############################################
## define the angle of the incident light
## if kx0 = 0 we only have to solve the eigenvalue problem once: time_compute = 1
## if kx0 != 0 we have to consider both the left and right eigenvector: time_compute = 2
kx0 = 0
time_compute = 1
# kx0 = 0.5*pi/d
# time_compute = 2
###############################################
par_gmsh_getdp = open('parameters_gmsh_getdp.dat', 'w')
par_gmsh_getdp.write('A           = %3.5e;\n'%(A          ))
par_gmsh_getdp.write('cel         = %3.5e;\n'%(cel        ))
par_gmsh_getdp.write('mu0         = %3.5e;\n'%(mu0        ))
par_gmsh_getdp.write('epsilon0    = %3.5e;\n'%(epsilon0   ))
par_gmsh_getdp.write('d           = %3.5e;\n'%(d          ))
par_gmsh_getdp.write('d_x         = %3.5e;\n'%(d_x        ))
par_gmsh_getdp.write('d_y         = %3.5e;\n'%(d_y        ))
par_gmsh_getdp.write('space2pml   = %3.5e;\n'%(space2pml  ))
par_gmsh_getdp.write('pmlsize     = %3.5e;\n'%(pmlsize    ))
par_gmsh_getdp.write('epsilon_oo_re  = %3.5e;\n'%(epsilon_oo.real))
par_gmsh_getdp.write('epsilon_oo_im  = %3.5e;\n'%(epsilon_oo.imag))
par_gmsh_getdp.write('omega_p     = %3.5e;\n'%(omega_p    ))
par_gmsh_getdp.write('gamma       = %3.5e;\n'%(gamma      ))
par_gmsh_getdp.write('eps_sup_re  = %3.5e;\n'%(eps_sup.real))
par_gmsh_getdp.write('eps_sup_im  = %3.5e;\n'%(eps_sup.imag))
par_gmsh_getdp.write('eps_sub_re  = %3.5e;\n'%(eps_sub.real))
par_gmsh_getdp.write('eps_sub_im  = %3.5e;\n'%(eps_sub.imag))
par_gmsh_getdp.write('npt_integ    = %d;\n'%(npt_integ       ))
par_gmsh_getdp.write('nb_slice     = %d;\n'%(nb_slice        ))
par_gmsh_getdp.write('scan_dist    = %3.5e;\n'%(scan_dist   ))
par_gmsh_getdp.write('ycut_sub_min = %3.5e;\n'%(ycut_sub_min))
par_gmsh_getdp.write('ycut_sub_max = %3.5e;\n'%(ycut_sub_max))
par_gmsh_getdp.write('ycut_sup_min = %3.5e;\n'%(ycut_sup_min))
par_gmsh_getdp.write('ycut_sup_max = %3.5e;\n'%(ycut_sup_max))
par_gmsh_getdp.write('kx0           = %3.5e;\n'%(kx0))
par_gmsh_getdp.write('expect_neig   = %3.5e;\n'%(expect_neig))
par_gmsh_getdp.write('eig_target_re = %3.5e;\n'%(eig_target_re))
par_gmsh_getdp.write('eig_target_im = %3.5e;\n'%(eig_target_im))
par_gmsh_getdp.write('angle = Pi/30;\n')
par_gmsh_getdp.write('fake_lambda  = 4e-2;;\n')
par_gmsh_getdp.close()
########################################################################################
########### MAIN ####################################################################
os.system('rm -rf field')
os.system('mkdir field')
# ##############################################
### create the mesh first 
os.system(str_gmsh_path+'gmsh '+mesh_geo+' -2 -o '+mesh_filename)

# ##### EIGENVALUE ####################################################################
slepc_options_st  = ' -st_type  sinvert -st_ksp_type preonly -st_pc_type lu -pc_factor_mat_solver_type mumps -petsc_prealloc 500 -v 0'
slepc_options_pep = ' -pep_max_it 400 -pep_target_magnitude'
slepc_options_rational = ' -nep_max_it 200 -nep_type nleigs -nep_rational -nep_target_magnitude -petsc_prealloc 50 -rg_interval_endpoints -10,0,0,50'
eigen_getdp1 = str_getdp_path+'getdp spectral1.pro -pre Modal1 -msh '+mesh_filename+' -cal -pos postop_Hz -slepc'                    
os.system(eigen_getdp1 + slepc_options_st+slepc_options_pep)
if time_compute == 2:
    eigen_getdp2 = str_getdp_path+'getdp spectral2.pro -pre Modal1 -msh '+mesh_filename+' -cal -pos postop_Hz -slepc'                    
    os.system(eigen_getdp2 + slepc_options_st+slepc_options_pep)

eigs = np.loadtxt('EigenValuesReal.txt',usecols=[5]) + 1j*np.loadtxt('EigenValuesImag.txt',usecols=[5])   
par_gmsh_getdp = open('parameters_gmsh_getdp.dat', 'a')
par_gmsh_getdp.write('neig = %3.5e;\n'%(len(eigs) )) 
# par_gmsh_getdp.write('neig = 1;\n') 
par_gmsh_getdp.close()

#### plot eigen mode #############
pl.figure()
vpr = np.loadtxt('./EigenValuesReal.txt',usecols=[5])
vpi = np.loadtxt('./EigenValuesImag.txt',usecols=[5])
pl.plot(vpr, vpi, 'yo')
for k in range(len(vpr)):
      pl.annotate(k, xy=(vpr[k],vpi[k]), xytext=None, xycoords='data', textcoords='data', arrowprops=None)
pl.title("Eigenvalues")
pl.xlabel("Re")
pl.ylabel("Im")
# # pl.figure()
# # vpr = np.loadtxt('./EigenValuesReal_ad.txt',usecols=[5])
# # vpi = np.loadtxt('./EigenValuesImag_ad.txt',usecols=[5])
# # pl.plot(vpr, vpi, 'yo')
# # for k in range(len(vpr)):
# #       pl.annotate(k, xy=(vpr[k],vpi[k]), xytext=None, xycoords='data', textcoords='data', arrowprops=None)
# # pl.title("Eigenvalues")
# # pl.xlabel("Re")
# # pl.ylabel("Im")
# # ###############################

##### norm ##################################
os.system('rm norm1.txt')
os.system('rm norm2.txt')
if time_compute == 2:
    norm_getdp = str_getdp_path+'getdp norm2.pro -pre Projection -msh  '+mesh_filename+' -cal -pos postop_norm -v 0'
else:  
    norm_getdp = str_getdp_path+'getdp norm.pro -pre Modal1 -msh  '+mesh_filename+' -res spectral1.res -pos postop_norm -v 0'
os.system(norm_getdp)

norm1 = np.loadtxt('norm1.txt',usecols=[5]) + 1j*np.loadtxt('norm1.txt',usecols=[6] )  # H*H all
norm2 = np.loadtxt('norm2.txt',usecols=[5]) + 1j*np.loadtxt('norm2.txt',usecols=[6] )  # dH*dH out
d1 = 2/cel**2 * eigs * norm1
d2 = ((2*eigs + 1j*gamma)*omega_p**2) / (epsilon_oo*(eigs**2+1j*eigs*gamma)-omega_p**2)**2 * norm2 
dd  = d1 + d2

#### LOOP FREQUENCY #################################################################
N_frequency = 15
N_kx = 2
R = np.zeros((N_kx,N_frequency))
T = np.zeros((N_kx,N_frequency))
R_ordre     = np.zeros((N_kx,N_frequency,2*N_d_order+1),dtype=complex)
T_ordre     = np.zeros((N_kx,N_frequency,2*N_d_order+1),dtype=complex)
Q_lamel_in  = np.zeros((N_kx,N_frequency))

## [32.5 38.2]
# lambda1 = 2*pi*cel/31
# tab_lambdas = np.linspace(lambda1,0.8,N_frequency)
# tab_lambdas = np.linspace(0.4,0.8,N_frequency)
tab_lambdas = np.linspace(0.5,0.8,N_frequency)

for i in range(N_frequency) :                   # loop over all the frequency
    lambda0 = tab_lambdas[i]
    theta = np.arcsin( kx0*lambda0*eps_sup.real/2/pi)
    omega0 = 2*pi*cel/lambda0   
    epsr_drude =  epsilon_oo - omega_p**2/(omega0**2+1j*omega0*gamma) 
    print (epsr_drude)
    par_gmsh_getdp = open('parameters_gmsh_getdp.dat', 'a')
    par_gmsh_getdp.write('lambda0 = %3.5e;\n'%(lambda0))
    par_gmsh_getdp.write('eps_lamel_in_re = %3.5e;\n'%(epsr_drude.real ))
    par_gmsh_getdp.write('eps_lamel_in_im = %3.5e;\n'%(epsr_drude.imag ))
    par_gmsh_getdp.close()

    #### COMPUTE PROJECTION PROBLEM #########################################################      
    ### compute Jn ######################################
    os.system('rm Jns.txt')
    if time_compute == 2:
        Jn_getdp = str_getdp_path+'getdp Jn.pro -pre Projection -msh '+mesh_filename+' -res spectral2.res -pos postop_Jn -v 0'
    else:
        Jn_getdp = str_getdp_path+'getdp Jn.pro -pre Projection -msh '+mesh_filename+' -res spectral1.res -pos postop_Jn -v 0'
    os.system(Jn_getdp)
    ### projection #####################################
    Jn  = np.loadtxt('Jns.txt',usecols=[5]) + 1j*np.loadtxt('Jns.txt',usecols=[6] )    
    Pns = Jn/((omega0-eigs)*dd)
    file_Pns = open('Pns.dat', 'w')
    # for k in range(len(eigs)):
    for k, Pn in enumerate(Pns):
        file_Pns.write('Pns_re_%d = %3.5e;\n'%(k,Pn.real))
        file_Pns.write('Pns_im_%d = %3.5e;\n'%(k,Pn.imag))
    file_Pns.close()
    project_getdp = str_getdp_path+'getdp projection.pro -pre Projection -msh  '+mesh_filename+' -res spectral1.res -pos  postop_modal -v 0'
    os.system(project_getdp)
    ### output file directory ###########################
    filelisttxt = [ f for f in os.listdir("./Views") if f.endswith(".txt")]
    filelistout = [ f for f in os.listdir("./Views") if f.endswith(".out")]
    for f in filelisttxt+filelistout:
        os.remove("./Views/"+f)
    if not os.path.exists("Views"):
        os.makedirs("Views")
    #####################################################
    projection_coeff_getdp = str_getdp_path+'getdp projection_coeff.pro -pre Projection -msh  '+mesh_filename+' -cal -pos postop_modal -v 0' 
    os.system(projection_coeff_getdp)
    [R0i,T0i,Ri,Ti,Qi_lamel_in] = diffraction_efficiency(str_gmsh_path,lambda0,theta,d,eps_sub,eps_sup,npt_integ,nb_slice,N_d_order)
    R[0,i] = Ri
    T[0,i] = Ti
    R_ordre[0,i,:] = R0i
    T_ordre[0,i,:] = T0i
    Q_lamel_in[0,i] = Qi_lamel_in

    #### COMPUTE DIRECT PROBLEM ####################################################################
    ### output file directory ###########################
    filelisttxt = [ f for f in os.listdir("./Views") if f.endswith(".txt")]
    filelistout = [ f for f in os.listdir("./Views") if f.endswith(".out")]
    for f in filelisttxt+filelistout:
        os.remove("./Views/"+f)
    if not os.path.exists("Views"):
        os.makedirs("Views")
    ####################################################
    direct_getdp = str_getdp_path+'getdp direct.pro -pre Scattering -msh  '+mesh_filename+' -cal -pos postop_scat -v 0' 
    os.system(direct_getdp)
    [R0i,T0i,Ri,Ti,Qi_lamel_in] = diffraction_efficiency(str_gmsh_path,lambda0,theta,d,eps_sub,eps_sup,npt_integ,nb_slice,N_d_order)
    R[1,i] = Ri
    T[1,i] = Ti
    R_ordre[1,i,:] = R0i
    T_ordre[1,i,:] = T0i
    Q_lamel_in[1,i] = Qi_lamel_in

# os.system(str_gmsh_path+'gmsh '+mesh_geo+' us.pos up.pos -merge compare.geo')
dumpglob2npz('output.npz',globals()) 
pl.show()

## END ##################################################################
