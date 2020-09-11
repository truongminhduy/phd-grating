################################################################################
#### PLOT REFECTION TRANMISSION ################################################
####################################################################!!!#########

import numpy as np
import scipy as sc
import scipy.io as scio
import pylab as plt
import pandas as pd
import sys
import seaborn as sns
# print(sys.version) 
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
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
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['figure.titlesize'] = 14


############ inpute color #################################################################
tab_colors = ['C0','C1','C2','C3','C4','C5','C6','C7','C8']
# tab_colors = ['b','r', 'g',  'y']

##### figure1 ##############################################################################
plt.figure(figsize=(12,12))
##### subplot #####################################################
gs = gridspec.GridSpec(5, 1)


##### subplot 1 #####################################
ax1 = plt.subplot(gs[2:,0])

### data for subplot 1a ##########
# res = np.load('output.npz')['saveddic'].item()
# res = np.load('save_k05/output_direct_pi20_a20_k05.npz')['saveddic'].item()
res = np.load('save_k05/output_direct_pi20_a11_k05_ref.npz')['saveddic'].item()
N_d_order = res['N_d_order']
T_ordre = res['T_ordre']
R_ordre = res['R_ordre']
Q_lamel_in   = res['Q_lamel_in']
N_kx = res['N_kx']
cel = res['cel']
tab_lambdas = res['tab_lambdas']
omega = 2*pi*cel/tab_lambdas
# compute total diffraction coefficients
Ttot = np.real(T_ordre.sum(axis=2))
Rtot = np.real(R_ordre.sum(axis=2))
Qtot = Q_lamel_in
Total = Rtot+Ttot+Qtot
### plot subplot 1a #############
i = 0
j = 0
ax1.plot(omega,Ttot[i,:],'-.',color=tab_colors[j],label='Transmission')
ax1.plot(omega,Rtot[i,:],'--',color=tab_colors[j],label='Reflection')
ax1.plot(omega,Qtot[i,:],'-',color=tab_colors[j],label='Absorption')
ax1.plot(omega,Total[i,:],':',color=tab_colors[j],label='T+R+A')

# ax1.plot(omega,Qtot[i,:],'-',color=tab_colors[j],label='Absorption')
# for mode in QNM[ee]:
# 	plt.axvline(mode)
plt.xlabel(r"$\operatorname{Re}(\omega_n)$ $(\times 10^{14} rad.s^{-1})$")
plt.ylabel("Diffraction efficiency")
plt.legend(loc =1)
# plt.xlim(omega[-1], omega[0])
plt.grid()
#####################################

### data for subplot 1b ##########
# res = np.load('output.npz')['saveddic'].item()
res = np.load('save_k05/output_pi20_a11_300_k05.npz')['saveddic'].item()
N_d_order = res['N_d_order']
T_ordre = res['T_ordre']
R_ordre = res['R_ordre']
Q_lamel_in   = res['Q_lamel_in']
tab_lambdas = res['tab_lambdas']
omega = 2*pi*cel/tab_lambdas
# compute total diffraction coefficients
Ttot = np.real(T_ordre.sum(axis=2))
Rtot = np.real(R_ordre.sum(axis=2))
Qtot = Q_lamel_in
Total = Rtot+Ttot+Qtot
### plot subplot 1b #############
i=0
j = 2
# lambdas = tab_lambdas*0.01
ax1.plot(omega,Ttot[i,:],'-.',color=tab_colors[j],label='Transmission')
ax1.plot(omega,Rtot[i,:],'--',color=tab_colors[j],label='Reflection')
ax1.plot(omega,Qtot[i,:],'-',color=tab_colors[j],label='Absorption')
ax1.plot(omega,Total[i,:],':',color=tab_colors[j],label='Trans+Refle')
################################

##### end subplot 1 #####################################


##### subplot 2 #########################################
ax2 = plt.subplot(gs[:2,0], sharex=ax1)

### data process
eigs = res['eigs']
vpr, vpi = eigs.real, eigs.imag 
lim1, lim2 = min(omega), max(omega)
# lim1, lim2 = min(vpr)*0.99, max(vpr)*1.01
plt.xlim(lim1, lim2)
plt.ylim(-0.05, 0)

##### plot eigen frequency ###############
scale_x = 1e-2
scale_y = 1e-2
ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_x))
ax2.xaxis.set_major_formatter(ticks_x)

ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_y))
ax2.yaxis.set_major_formatter(ticks_y)

plt.plot(vpr, vpi, 'C1o',markersize=4 ,mfc='none')
# for k in range(len(vpr)):
#       plt.annotate(k, xy=(vpr[k],vpi[k]), xytext=None, xycoords='data', textcoords='data', arrowprops=None)
# plt.title("Eigenvalues")
# plt.xlabel("Re")
plt.ylabel(r"$\operatorname{Im}(\omega_n)$ $(\times 10^{14} rad.s^{-1})$")
ax2.xaxis.grid()

### set axis for eigen frequency ###############
ax3 = ax2.twiny()
freq = np.linspace(lim1, lim2, 10)
fake = np.linspace(0, 0, 10)
lamb = np.round(2*pi*cel/freq*0.01,2)
### plot invisible points ############
ax3.plot(freq, fake, linestyle = "None", marker = "o", color = 'w', markersize = 0.01)
### change label of axis
plt.xticks(freq, lamb)
ax3.set_xlim(ax2.get_xlim())
ax3.set_xlabel(r"$\lambda (\mu m)$")
##### end subplot 1 #####################################


plt.show()



