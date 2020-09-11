###############################################################################
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
# ax1 = plt.figure(figsize=(12,12))
# plt.figure()
ax1 = plt.subplot(111)
##### subplot #####################################################
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
# ax1.plot(omega,Ttot[i,:],'-.',color=tab_colors[j],label='Transmission')
# ax1.plot(omega,Rtot[i,:],'--',color=tab_colors[j],label='Reflection')
# ax1.plot(omega,Total[i,:],':',color=tab_colors[j],label='T+R+A')
ax1.plot(omega,Qtot[i,:],'-',color=tab_colors[j],label='direct')

plt.xlabel(r"$\operatorname{Re}(\omega_n)$ $(\times 10^{14} rad.s^{-1})$")
plt.ylabel('Absorption efficiency')
plt.legend(loc =1)

plt.grid()


########################
res = np.load('save_k05/output_pi20_a11_1_k05.npz')['saveddic'].item()
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
j = 3
# ax1.plot(omega,Ttot[i,:],'-.',color=tab_colors[j],label='Transmission')
# ax1.plot(omega,Rtot[i,:],'--',color=tab_colors[j],label='Reflection')
# ax1.plot(omega,Total[i,:],':',color=tab_colors[j],label='T+R+A')
ax1.plot(omega,Qtot[i,:],'-',color=tab_colors[j],label='1')
plt.legend(loc =1)
######################
res = np.load('save_k05/output_pi20_a11_50_k05.npz')['saveddic'].item()
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
j = 4
# ax1.plot(omega,Ttot[i,:],'-.',color=tab_colors[j],label='Transmission')
# ax1.plot(omega,Rtot[i,:],'--',color=tab_colors[j],label='Reflection')
# ax1.plot(omega,Total[i,:],':',color=tab_colors[j],label='T+R+A')
ax1.plot(omega,Qtot[i,:],'-',color=tab_colors[j],label='50')
plt.legend(loc =1)


#####################
res = np.load('save_k05/output_pi20_a11_180_k05.npz')['saveddic'].item()
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
j = 5
# ax1.plot(omega,Ttot[i,:],'-.',color=tab_colors[j],label='Transmission')
# ax1.plot(omega,Rtot[i,:],'--',color=tab_colors[j],label='Reflection')
# ax1.plot(omega,Total[i,:],':',color=tab_colors[j],label='T+R+A')
ax1.plot(omega,Qtot[i,:],'-',color=tab_colors[j],label='180')
plt.legend(loc =1)

#####################
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
j = 8
# ax1.plot(omega,Ttot[i,:],'-.',color=tab_colors[j],label='Transmission')
# ax1.plot(omega,Rtot[i,:],'--',color=tab_colors[j],label='Reflection')
# ax1.plot(omega,Total[i,:],':',color=tab_colors[j],label='T+R+A')
ax1.plot(omega,Qtot[i,:],'-',color=tab_colors[j],label='300')
plt.legend(loc =1)

scale_x = 1e-2
ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_x))
ax1.xaxis.set_major_formatter(ticks_x)
plt.show()

