#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
from scipy.constants import golden 
import os
import numpy as np
import pandas as pd

# Plotting parameters, to make it nicer
matplotlib.rcParams['figure.figsize'] = (6.5 * golden,6.5)
matplotlib.rcParams['figure.subplot.wspace'] = 0.5
matplotlib.rcParams['figure.subplot.left'] = 0.1
matplotlib.rcParams['figure.subplot.hspace'] = 0.3
matplotlib.rcParams['figure.subplot.bottom'] = 0.2
matplotlib.rcParams['figure.subplot.right'] = 1
matplotlib.rcParams['figure.subplot.top'] = 0.85
matplotlib.rcParams['axes.labelsize'] = 22
matplotlib.rcParams['axes.titlesize'] = 22
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.size'] = 22

## Mechanical functions

# Compute Biot modulus
def Biot_modulus(n, K, a, K_f):
    Biot_modulus = n / K_f + ((1.0 - a) * (a - n)) / K
    return Biot_modulus

# Compute compressibility
def Confined_compressibility(K, G):
    Confined_compressibility = 1.0 / (K + (4.0 / 3.0) * G)
    return Confined_compressibility

# Compute consolidation coefficient
def Consolidation_coeff(k, a, y, Biot_modulus, Confined_compressibility):
    Consolidation_coeff = k / (y * (Biot_modulus + a**2.0 * Confined_compressibility)) 
    return Consolidation_coeff

# Compute initial condition (given by compression of the pore structure)
def Initial_condition(a, Biot_modulus, Confined_compressibility, q):
    Initial_condition = ((a * Confined_compressibility) / (Biot_modulus + a ** 2.0 * Confined_compressibility)) * q 
    return Initial_condition

# Terzaghi's problem solution
def Terzaghis_problem(c_v, t, h, n_serie = None):
    
    ##--------------------------
    # Terzaghis_problem function:
    # Function to compute classical Terzaghi's problem of consolidation of soils
        # input
            # c_v: Consolidation coefficient (SI units) Float
            # t: Time (SI units) List
            # h: Depth (SI units) Float
            # n_serie: Optional value with the number of eigenvalues computed (the more the better) Integer
        
        # return
            # pressure: Array ([101, length of t])
            # Normalized lenght: lenght discretization / Total lenght
    ##--------------------------

    # Make time an array so we can looped over it later
    t = np.array(t)
    
    # Points where the solution is computed (change this parameter if you want more points un your solution)
    Z = np.linspace(0, h, 101) 
    
    # Empty array to store pressure solution. Rows = length scale, Columns = time 
    pressure = np.zeros((len(Z),len(t))) 
    
    n_serie = n_serie or 100
    
    # Time loop. Time, element form t list. j is used to append results in pressure array
    for j, time in enumerate(t):
        
        # Lenght loop. z is depth discretization, i is used to append results in pressure array
        for i, z in enumerate(Z):
            # Empty array to store sum of the loop below
            D = np.zeros(n_serie - 1) 

            # Compute sum from the inverse Fourier transfor solution. A, B, C are used for simplification.
            # solution is storaged in D, then the sum of the array is obtained.
            for k in range(1, n_serie): 
                A = ((-1.0) ** (k - 1.0)) / (2.0 * k - 1.0)
                B = np.cos((2.0 * k - 1.0) * (np.pi / 2.0) * (z / h))
                C = np.exp( - ((2.0 * k - 1.0) ** 2.0) * ((np.pi ** 2.0) / 4.0) * (( c_v * time) / (h ** 2.0)))
                D[k - 1] = (A * B * C)
            
            # Get solution
            pressure[i, j] = ((4.0 / np.pi) * D.sum())

            
    return pressure, -Z[::-1]/h

# System's characteristics
rho = 1000 # Fluid's density
nu = 1e-3 # Fluid's viscocity
g = 10 # Gravity's acceleration
Kf = 2.2E9 # Fluid's bulk modulus
a = 0.9 # Biot's coefficient

# Analytical solution

# Hydraulic parameters of the soil
k = 1E-4 # Permeability
n = 0.2 # Porosity

# Mechanical parameters of the soil
L = 8.4E7 # Bulk modulus
G = 6.25E7 # Shear modulus

# Water parameters
y = 1E4 # Volumetric weigth of fluid

#System's characteristics
q = 1000 # Load
h = 100 # Depth of the system

# Compute poroelasticity constants 
S = Biot_modulus(n, L, a, Kf) 
m_v = Confined_compressibility(L, G) 
c_v = Consolidation_coeff(k, a, y, S, m_v)
p_0 = Initial_condition(a, S, m_v, q)

# Compute solution for given times in s
t = [1, 10, 50, 200, 1000]

# Call function
p, z = Terzaghis_problem(c_v, t, h)

# Plot
fig, ax = plt.subplots()

for i, time in enumerate(t):
    ax.plot(p[:,i], np.linspace(0,1,len(p)),  linewidth=2, label='t = ' + str(time) + 's')

# Numerical solution
df = pd.read_csv('gold/TerzaghiImportData_full.csv')
idx = [0, 1, 2, 3, 4] # Indices where the solution is found
point_val = np.linspace(0, 1, 11) # Space array
for i in idx:
    plt.plot(df.iloc[i][1:] / p_0, point_val, 'k.', markersize=5, label='RHEA\'s solution' if i == 0 else '')

ax.legend(fontsize=20,ncol=2, handleheight=1.5, labelspacing=.1)
ax.set_xlabel('$p_f/p_0$')
ax.set_ylabel('z/L');

plt.savefig('Terzaghis_problem.png', bbox_inches='tight')
