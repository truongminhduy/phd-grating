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

### MAIN ################################################################
res = np.load('save_k00_1/output_pi30_2a25_720.npz')['saveddic'].item()

epsilon_oo = res['epsilon_oo']
omega_p    = res['omega_p']
gamma      = res['gamma']
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
Total1 = Rtot+Ttot+Qtot

res = np.load('save_k00_1/output_pi20_2a25_750.npz')['saveddic'].item()
T_ordre    = res['T_ordre']
R_ordre    = res['R_ordre']
Q_lamel_in = res['Q_lamel_in']
tab_lambdas = res['tab_lambdas']
Ttot = np.real(T_ordre.sum(axis=2))
Rtot = np.real(R_ordre.sum(axis=2))
Qtot = Q_lamel_in
Total2 = Rtot+Ttot+Qtot

res = np.load('save_k00_1/output_pi15_2a25_750.npz')['saveddic'].item()
T_ordre    = res['T_ordre']
R_ordre    = res['R_ordre']
Q_lamel_in = res['Q_lamel_in']
tab_lambdas = res['tab_lambdas']
Ttot = np.real(T_ordre.sum(axis=2))
Rtot = np.real(R_ordre.sum(axis=2))
Qtot = Q_lamel_in
Total3 = Rtot+Ttot+Qtot

# res = np.load('sav_k00_1/output_pi30_2a30_1200.npz')['saveddic'].item()
# T_ordre    = res['T_ordre']
# R_ordre    = res['R_ordre']
# Q_lamel_in = res['Q_lamel_in']
# tab_lambdas = res['tab_lambdas']
# Ttot = np.real(T_ordre.sum(axis=2))
# Rtot = np.real(R_ordre.sum(axis=2))
# Qtot = Q_lamel_in
# Total0 = Rtot+Ttot+Qtot


plt.figure(figsize=(7,5.5))
# plt.figure()
ax = plt.subplot(111);
# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width*0.85, box.height])
err1 = abs(Total1[0,:]-1)
err2 = abs(Total2[0,:]-1)
err3 = abs(Total3[0,:]-1)
# ax.plot(omegas,Total1[1,:],'C0:',label=r'$direct$')
ax.plot(omegas,err1,'C0-', linewidth=2,label=r'$\phi=\pi/30$')
ax.plot(omegas,err2,'C1-.',linewidth=2,label=r'$\phi=\pi/20$')
ax.plot(omegas,err3,'C2--',linewidth=2,label=r'$\phi=\pi/15$')
plt.xlabel(r'$\operatorname{Re}(\omega)$')
plt.ylabel("Total energy - Absolute error")
plt.legend()

plt.grid()
# plt.savefig('energy.eps', format='eps', dpi=1000)
# plt.savefig('energy.png', format='png', dpi=1000)
plt.show()
