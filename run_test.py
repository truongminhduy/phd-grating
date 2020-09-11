import numpy as np
import scipy as sc
import scipy.io as scio
import pylab as pl
import os
import sys
import scipy.signal as signal
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

np.set_printoptions(precision=3)
pi = np.pi
#### MAIN 1 #######################################################################
### define inpute files ###########################################################
# path = '../../../Downloads/gmsh-4.3.0-Linux64/bin/'
path = ''
# os.system(path+'gmsh t1.geo')
os.system(path+'gmsh mesh_duy_simple.geo ')
# os.system(path+'gmsh mesh_duy_simple.geo us.pos up.pos -merge compare.geo')
