import numpy as np
import scipy as sc
import scipy.io as scio
import pylab as plt
import pandas as pd
import sys
import seaborn as sns
# print(sys.version) 

import os
import sys
import matplotlib.pyplot as plt
plt.rc('font', size=13,)
pi = np.pi
plt.style.use('seaborn-muted')
plt.rc('text', usetex=True)
plt.rc('font', **{'family' : "serif"})
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(params)
#plt.rcParams['text.latex.preamble'] = [r'\boldmath']
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['figure.titlesize'] = 14


# res = np.load('output.npz')['saveddic'].item()
# res = np.load('save_k00_1/output_pi30_2a30_1200.npz')['saveddic'].item()
res = np.load('save_k00_1/output_pi30_2a25_720.npz')['saveddic'].item()
epsilon_oo = res['epsilon_oo']
omega_p    = res['omega_p']
gamma      = res['gamma']
eigs       = res['eigs']
cel        = res['cel']

print (2*pi*cel/0.4)

pole_root = np.array([0, -1j*gamma])
plt.figure()
plt.plot(eigs.real, eigs.imag, 'C0o',markersize=3 ,mfc='none',label='Eigenmode10')
plt.xlabel(r"$\operatorname{Re}(\omega)$") 
plt.ylabel(r"$\operatorname{Im}(\omega)$")


### POLE PROBLEM ##################
plt.plot(pole_root.real,pole_root.imag,'C2x',markersize=6,label='Pole')

# ### CORNER PROBLEM ###############
# k = np.linspace(-3,-1/3,20)
# omega_k = -1j*gamma/2 + np.sqrt(gamma**2/4+omega_p**2/(1-k))
# plt.plot(omega_k.real, omega_k.imag, 'C2-',markersize=5 ,mfc='none',label='asd')

# ### ESPI=0 PROBLEM ##################
# omega_0 = -1j*gamma/2 + np.sqrt(gamma**2/4+omega_p**2/epsilon_oo)
# plt.plot(omega_0.real,omega_0.imag,'C2x',label='E=0')


plt.legend()
plt.show()