import numpy as np
import scipy as sc
import os
pi = np.pi

def diffraction_efficiency(str_gmsh_path,lambda0,theta,d,eps_sub,eps_sup,npt_integ,nb_slice,N_d_order):
    decalage  = 0
    No_ordre  = np.linspace(-N_d_order+decalage, N_d_order+decalage, 2*N_d_order+1)
    # x_slice   = np.arange(-d/2, d/2,d/npt_integ)
    npt_integ = int(npt_integ)
    x_slice   = sc.linspace(-d/2, d/2, npt_integ)
    k_sub     = 2*pi*sc.sqrt(eps_sub)/lambda0
    k_sup     = 2*pi*sc.sqrt(eps_sup)/lambda0
    alpha_sup = k_sup*sc.sin(theta)                
    beta_sup  = sc.sqrt(k_sup**2-alpha_sup**2)
    beta_sub  = sc.sqrt(k_sub**2-alpha_sup**2)
    s_t		  = sc.zeros((1,(2*N_d_order+1)),complex)[0,:]
    s_r		  = sc.zeros((1,(2*N_d_order+1)),complex)[0,:]
    Aeff_t    = sc.zeros((nb_slice,2*N_d_order+1),complex)
    Aeff_r	  = sc.zeros((nb_slice,2*N_d_order+1),complex)
    Hz_diff_t = np.loadtxt('./Views/sub_field_cuts.out',usecols=[8])+1j*np.loadtxt('./Views/sub_field_cuts.out',usecols=[9])
    Hz_diff_r = np.loadtxt('./Views/sup_field_cuts.out',usecols=[8])+1j*np.loadtxt('./Views/sup_field_cuts.out',usecols=[9])
    Hz_diff_t = np.transpose(Hz_diff_t.reshape(npt_integ,nb_slice,order="F"))
    Hz_diff_r = np.transpose(Hz_diff_r.reshape(npt_integ,nb_slice,order="F"))
    for m1 in range(0,nb_slice):
        slice_t   = Hz_diff_t[m1,:]         # solution u
        slice_r   = Hz_diff_r[m1,:]
        alphat_t  = alpha_sup + 2*pi/(d)*No_ordre
        alphat_r  = alpha_sup + 2*pi/(d)*No_ordre
        betat_sup = sc.sqrt(k_sup**2-alphat_r**2)
        betat_sub = sc.sqrt(k_sub**2-alphat_t**2)
        for k in range(0,2*N_d_order+1):
            expalpha_t = sc.exp(-1j*alphat_t[k]*x_slice);
            expalpha_r = sc.exp(-1j*alphat_r[k]*x_slice);
            s_t[k] = sc.trapz(slice_t*expalpha_t,x=x_slice)/d
            s_r[k] = sc.trapz(slice_r*expalpha_r,x=x_slice)/d
        Aeff_t[m1,:] = (np.abs(s_t))**2*betat_sub/beta_sup/eps_sub          
        Aeff_r[m1,:] = (np.abs(s_r))**2*betat_sup/beta_sup
    Rordre = np.mean(Aeff_r,axis=0)
    Tordre = np.mean(Aeff_t,axis=0)
    R = np.mean(np.sum(Aeff_r.real,axis=1))
    T = np.mean(np.sum(Aeff_t.real,axis=1))
    Q_lamel_in = np.loadtxt('./Views/temp-Q_lamel_in.txt')[1]
    # Q_lamel_in = 0
    # print ('\n******************')
    # print ('* ENERGY BALANCE *'  )
    # print ('******************'  )
    # print ('R             = ', "%0.6f"%R ,'     (standard dev slice2slice=',sc.std(np.sum(Aeff_r.real,axis=1)),')')
    # print ('T             = ', "%0.6f"%T ,'     (standard dev slice2slice=',sc.std(np.sum(Aeff_t.real,axis=1)),')')
    # print ('Q_lamel_in    = ', "%0.6f"%Q_lamel_in)
    # # print 'Q_layer_acc   = ', "%0.6f"%Q_layer_acc
    # print ('------------------------')
    # print ('TOTAL        = ', "%0.6f"%(T+R+Q_lamel_in))
    return [Rordre,Tordre,R,T,Q_lamel_in]
