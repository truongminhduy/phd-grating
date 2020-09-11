################################################################################
#### PLOT REFECTION TRANMISSION ################################################
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
plt.rcParams['font.size'] = 20
plt.rcParams['axes.labelsize'] = 21
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['figure.titlesize'] = 20
plt.rcParams['lines.linewidth'] = 2
### MAIN ################################################################
# res = np.load('output.npz')['saveddic'].item()
# res = np.load('save_k00_1/output_pi15_2a20_800.npz')['saveddic'].item()
# res = np.load('save_k00_1/output_pi20_2a20_800.npz')['saveddic'].item()
# res = np.load('save_k00_1/output_pi15_2a20_800.npz')['saveddic'].item()
# res = np.load('save_k00_1/output_pi20_2a25_750.npz')['saveddic'].item()
res = np.load('save_k00_1/output_pi30_2a25_720.npz')['saveddic'].item()
# res = np.load('save_k00_1/output_pi30_2a30_1200.npz')['saveddic'].item()

epsilon_oo = res['epsilon_oo']
omega_p    = res['omega_p']
gamma      = res['gamma']
# eigs       = res['eigs']
cel        = res['cel']
N_d_order  = res['N_d_order']
T_ordre    = res['T_ordre']
R_ordre    = res['R_ordre']
Q_lamel_in = res['Q_lamel_in']
tab_lambdas = res['tab_lambdas']
omegas = 2*pi*cel/tab_lambdas

Ttot = np.real(T_ordre.sum(axis=2))
Rtot = np.real(R_ordre.sum(axis=2))
Qtot = Q_lamel_in
Total = Rtot+Ttot+Qtot

plt.figure(figsize=(12,8))
ax = plt.subplot(111);
# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width*0.85, box.height])
i = 1
tab_colors = 'C0'
ax.plot(omegas,Ttot[i,:],'-.',color=tab_colors,label='Transmission')
ax.plot(omegas,Rtot[i,:],'--',color=tab_colors,label='Refection')
ax.plot(omegas,Qtot[i,:],'-',color=tab_colors,label='Absorption')
ax.plot(omegas,Total[i,:],':',color=tab_colors,label='R+T+A')
plt.xlabel(r'$\operatorname{Re}(\omega)$')
plt.ylabel("Diffraction efficiency")
plt.legend(loc=6)

# res = np.load('output.npz')['save_k00_1ddic'].item()
# res = np.load('save_k00_1/output_pi30_2a25_720.npz')['saveddic'].item()
cel        = res['cel']
N_d_order  = res['N_d_order']
T_ordre    = res['T_ordre']
R_ordre    = res['R_ordre']
Q_lamel_in = res['Q_lamel_in']
tab_lambdas = res['tab_lambdas']
omegas = 2*pi*cel/tab_lambdas
Ttot = np.real(T_ordre.sum(axis=2))
Rtot = np.real(R_ordre.sum(axis=2))
Qtot = Q_lamel_in
Total = Rtot+Ttot+Qtot
i = 0
tab_colors = 'C2'
ax.plot(omegas,Ttot[i,:],'-.',color=tab_colors,label='T ')
ax.plot(omegas,Rtot[i,:],'--',color=tab_colors,label='R')
ax.plot(omegas,Qtot[i,:],'-',color=tab_colors,label='Q')
ax.plot(omegas,Total[i,:],':',color=tab_colors,label='R+T+Q')


plt.grid()

# plt.savefig('deff.eps', format='eps', dpi=1000)
# plt.savefig('deff.png', format='png', dpi=1000)
plt.show()
