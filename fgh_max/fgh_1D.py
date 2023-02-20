#
# Numerical solution of the 1D Schroedinger equation in an arbitrary bound potential
# Fourier Grid Hamiltonian method:
#
# [1] Marston, C. C.; Balint-Kurti, G. G.
#     The Fourier Grid Hamiltonian Method for Bound State Eigenvalues and Eigenfunctions.
#     J. Chem. Phys. 1989, 91, 3571–3576
#     https://doi.org/10.1063/1.456888
#
# [2] Balint-Kurti, G. G.; Ward, C. L.; Marston, C. C.
#     Two Computer Programs for Solving the Schrödinger Equation for Bound-State Eigenvalues
#     and Eigenfunctions Using the Fourier Grid Hamiltonian Method.
#     Computer Physics Communications 1991, 67, 285–292
#     https://doi.org/10.1016/0010-4655(91)90023-e
#
# Original code by Maxim Secor, Yale University, 2022
#

import math
import numpy as np
import os
import sys
from numba import jit

#from scipy.interpolate import interp1d
#import matplotlib.pyplot as plt

# mass - mass of the particle in atomic units (electron mass units, proton's mass is 1836.15)
# potential - array of the portential energy values (in Hartrees) at the equally spaced grid points
# grid - array of EQUALLY SPACED grid points (in Bohrs)
#       (the more the better, for a proton ~128 points within ~2A is usually enough for acceptable accuracy.
#        The results should be convergent with respect to the number of grid points)
#
# returns: eigenvalues (energy levels) and normalized eigenvectors (wavefunction values at the grid points)
# https://numpy.org/doc/stable/reference/generated/numpy.linalg.eigh.html
#

# the following makes the code much faster
@jit(nopython=True)

def fgh_1D(grid,potential,mass):

   nx = len(grid)
   dx = grid[1]-grid[0]
   k = np.pi/dx

   vmat = np.zeros((nx,nx))
   tmat = np.zeros((nx,nx))
   hmat = np.zeros((nx,nx))

   for i in range(nx):
       for j in range(nx):

           if i == j:
               vmat[i,j] = potential[j]
               tmat[i,j] = (k**2)/3
           else:
               dji = j-i
               vmat[i,j] = 0
               tmat[i,j] = (2*k**2)/(np.pi**2)*(((-1)**dji)/(dji**2))

           hmat[i,j] = (1/(2*mass))*tmat[i,j] + vmat[i,j]

   hmat_soln = np.linalg.eigh(hmat)
   return hmat_soln

#%%

