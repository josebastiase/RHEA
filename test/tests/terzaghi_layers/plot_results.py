#!/usr/bin/env python3
import numpy as np
from scipy.optimize import fsolve
from pychebfun import *
from functools import partial
from scipy import special
import scipy.integrate as integrate
from scipy.special import erf
import matplotlib
import matplotlib.pyplot as plt
from scipy.constants import golden 
import os
import pandas as pd
import matplotlib.ticker as mticker
fh = mticker.ScalarFormatter(useOffset=False, useMathText=True)
gh = lambda x,pos : "${}$".format(fh._formatSciNotation('%1.10e' % x))
fmt = mticker.FuncFormatter(gh)

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

# Terzaghi's problem solution multilayers
def Terzaghis_problem_multi(c_v, t, h, m, l, a, b, c, init, n_serie = None, space_sol = None):
    
    ##--------------------------
    # Terzaghis_problem_multi function:
    # Function to obtain the multi layers solution of Terzaghi's problem. The solution was obtained by 
    # Hickson, R. I., Barry, S. I., & Mercer, G. N. (2009). Critical times in multilayer diffusion. 
    # Part 1: Exact solutions. International Journal of Heat and Mass Transfer, 52(25-26), 5776-5783.
        # input
            # c_v: Consolidation coefficient (SI units) Float
            # t: Time (SI units) List
            # h: Depth of each layer (SI units) List
            # m: Number of layers Integer
            # l: Lenght of the layers (SI units) List
            # a, b, c : Muxed boundary conditions. Further information can be found in the paper above
            # n_serie: Optional value with the number of eigenvalues computed (the more the better) Integer
            # space_sol: Space discretization
        
        # return
            # pressure: Array ([space_sol, length of t])
            # Normalized lenght: lenght discretization
    ##--------------------------
    
    # Time properties
    t = np.array(t)
    init = np.array(init)
    
    # System properties
    l = np.array(l)
    z = np.zeros(m)
    
    # z stores the SS solution at the right boundaries 
    for i in range(m):
        if i == 0:
            z[0] = l[0]
        else:
            z[i] = z[i - 1] + l[i]
    
    # Diffusion properties. Notation D is used from diffusion coefficient
    D = np.array(c_v)
    D_av = sum(l / D)
    d = np.sqrt(D)
    
    # Boundary conditions
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)
    
    eigen_val = n_serie or 200
    space_sol = space_sol or m * 100
    
    # Compute steady state solution
    q = np.zeros(m)
    h = np.zeros(m)
    l_D = np.zeros(m)

    q[0] = ((a[0] * c[1] - a[1] * c[0]) * D[-1]) / \
    (a[0] * b[1] * D[0] - a[1] * b[0] * D[-1] + (a[0] * a[1] * D[0] * D[-1]) * D_av)

    h[0] = (c[0] - b[0] * q[0]) / (a[1] + 1E-8) # Compute h_0

    # Compute vectorized h and q
    for i in range(1, m):
        q[i] = D[0] * q[0] / D[i] 
        l_D[i] = l[i - 1] / D[i - 1]
        h[i] = h[0] + D[0] * q[0] * sum(l_D)
        
    # SS solution g    
    def steady_solution(i, x_i, x):
        g = init[i] - (q[i] * (x - x_i) + h[i])

        return g
    
    # SS solution w
    def w_solution(i, x_i, x):
        w = (q[i] * (x - x_i) + h[i])

        return w
    
   
    # Eigenfunction solution length-depend. mu_ is a eigenvalue and mu is the eigenvector
    def eigen_length(i, d, mu_, K1i, K2i, x_i, x): 
        params = K1i * np.sin((mu_ / d[i]) * (x - x_i)) + \
        K2i * np.cos((mu_ / d[i]) * (x - x_i)) 
    
        return params

    
    # Eigenfunction solution time-depend. K1i1 and K2i1 include other functions that are called later on
    def eigen_time(t, mu_):
        params = np.exp(- mu_ ** 2.0 * t)
        
        return params
    
    
    # Define K functions K_11 is always 1
    def K11(mu_):
    
        return 1


    def K21(b, a, d, mu_):
        params = - b[0] * mu_ / (a[0] * d[0])

        return params


    def K1i1(i, d, l, K1i, K2i, mu_):
        params = (d[i - 1] / d[i]) * (K1i(mu_) * np.cos(mu_ * l[i - 1] /\
        d[i - 1]) - K2i(mu_) * np.sin(mu_ * l[i - 1] / d[i - 1]))

        return params


    def K2i1(i, d, l, K1i, K2i, mu_):
        params = K1i(mu_) * (np.sin(mu_ * l[i - 1] / d[i - 1])) + K2i(mu_) * \
        (np.cos(mu_ * l[i - 1] / d[i - 1]))

        return params

    
    # Save K functions in an array. The purpose of this is to find K(mu_)
    K = np.empty((2, m), dtype=type(K2i1)) # Empty array to save the functions
    
    K[0, 0] = K11 # Save first function K_11
    K[1, 0] = partial(K21, b, a, d) # Partial is a function wrapper

    # After the first two are saved the others can be looped using partial
    for i in range(1, m):
        K[0, i] = partial(K1i1, i, d, l)  
        K[1, i] = partial(K2i1, i, d, l)  

    for n in range(1, m):
        K[0, n] = partial(K[0, n], K[0, n - 1], K[1, n - 1])
        K[1, n] = partial(K[1, n], K[0, n - 1], K[1, n - 1])

        
    
    # Function to wrap K_n to get K_n(mu_). Now eigenvalues are obtained
    def eigen_wrapper(K, d, l):

        def eigen(mu_):
            K1n = K[0, K.shape[1] - 1](mu_) # Get k1n to put it in the params function
            K2n = K[1, K.shape[1] - 1](mu_) # Get k2n to put it in the params function
            params = K1n * (a[1] * np.sin(mu_ * l[-1] / d[-1]) + \
            (mu_ * b[1] / d[-1]) * np.cos(mu_ * l[-1] / d[-1])) + \
            K2n * ((-mu_ * b[1] / d[-1]) * np.sin(mu_ * l[-1] / d[-1]) + \
            a[1] * np.cos(mu_ * l[-1] / d[-1]))
            return params

        return eigen
    
    # Get eigenvalues. Chebfunctions are used - Faster and better
    eigen_fun = chebfun(eigen_wrapper(K, d, l), [0, n_serie])
    mu = eigen_fun.roots() # eigenvector
    
    # Compute Sturm–Liouville and store the top solution in SLt and the bottom in SLb
    SLt = np.zeros((len(mu), m))
    SLb = np.zeros((len(mu), m))
    
    # slow part of the code. Get the sum of the integrals 
    for i in range(m):
        if i == 0:
            x_i = h[0]
            
            g_fun = chebfun(partial(steady_solution, i, x_i), [x_i, z[i]])
    
            for j, mu_ in enumerate(mu):
                K1i = K[0, i](mu_)
                K2i = K[1, i](mu_)
                eigen_length_fun = chebfun(partial(eigen_length, i, d, mu_, K1i, K2i, x_i), [x_i, z[i]])

                SLb[j, i] = eigen_length_fun.dot(eigen_length_fun)
                SLt[j, i] = g_fun.dot(eigen_length_fun)
                
        else:
            x_i = z[i - 1]
            g_fun = chebfun(partial(steady_solution, i, x_i), [x_i, z[i]])
            
            for j, mu_ in enumerate(mu):
                K1i = K[0, i](mu_)
                K2i = K[1, i](mu_)
                eigen_length_fun = chebfun(partial(eigen_length, i, d, mu_, K1i, K2i, x_i), [x_i, z[i]])

                SLb[j, i] = eigen_length_fun.dot(eigen_length_fun)
                SLt[j, i] = g_fun.dot(eigen_length_fun)

                

    SLt=SLt.sum(axis=1) 
    SLb=SLb.sum(axis=1)

    x = np.linspace(h[0], z[-1], int(space_sol))
    x_layer = int(space_sol / m)
    pressure = np.zeros((len(x),len(t)))

    #Time loop 
    for v, time in enumerate(t):

        #Layer loop
        for i in range(m):
            if i==0:
                x_i = h[0]
                x_disp = np.linspace(x_i, z[i], x_layer)

                #Lenght loop
                for j, r in enumerate(x_disp):
                    P = np.zeros_like(mu)
                    w = w_solution(i, x_i, r)

                    #Eigenvalues
                    K1i = K[0, i](mu)
                    K2i = K[1, i](mu)
                    time_fun = eigen_time(time, mu)
                    xi = eigen_length(i, d,mu, K1i, K2i, x_i, r)
                    P = (time_fun * xi * SLt / SLb)

                    pressure[j, v] = w + P.sum()

            else:
                x_i = z[i - 1]
                x_disp = np.linspace(x_i, z[i], x_layer)

                #Lenght loop
                for j, r in enumerate(x_disp):
                    P = np.zeros_like(mu)
                    w = w_solution(i, x_i, r)

                    #Eigenvalues 
                    K1i = K[0, i](mu)
                    K2i = K[1, i](mu)
                    time_fun = eigen_time(time, mu)
                    xi = eigen_length(i, d, mu, K1i, K2i, x_i, r)
                    P = (time_fun * xi * SLt / SLb)

                    pressure[(i * 100) + j, v] = w + P.sum()
                    
                    
    return pressure, x


# System's characteristics
rho = 1000 # Fluid's density
nu = 1e-3 # Fluid's viscocity
g = 10 # Gravity's acceleration
Kf = 2.2E9 # Fluid's bulk modulus
a = 0.9 # Biot's coefficient

# Analytical solution
m = 10

# Hydraulic parameters of the soil
k = np.linspace(1E-4, 1E-4, m) # Permeability
k[0::2] = 1E-8

n = np.linspace(.2, .2, m) # Porosity


# Mechanical parameters of the soil
L = np.linspace(8.4E7,8.4E7, m) # Bulk modulus of the soil

G = np.linspace(6.25E7, 6.25E7, m) # Shear modulus of the soil

# Water parameters
y = 1.0E4 # Volumetric weigth of fluid

q = 1000 # Load

# Compute poroelasticity constants
S = Biot_modulus(n, L, a, Kf)
m_v = Confined_compressibility(L, G)
c_v = Consolidation_coeff(k, a, y, S, m_v) # Consolidation coeffcient
p_0 = Initial_condition(a, S, m_v, q) # Initial conditions in each layer

t = [1E4, 5E6, 2E6, 7E6] # Time where the solution is computed
h = [0, 100] # Top and bottom of your model # Number of layers
l = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10] # Length of the layers
ab = [1.0 , 0.0] # Mixed boundary conditions, Dirichet
b = [0.0, 1.0] # Mixed boundary conditions, Neumann
c = [0.0, 0.0] # Mixed boundary conditions, value

p, z = Terzaghis_problem_multi(c_v, t, h, m, l, ab, b, c, p_0, n_serie = .05, space_sol = None)

# Plot

fig, ax = plt.subplots()
for i, time in enumerate(t):
    ax.plot(p[:,i]/p_0[0], z[::-1]/100, linewidth=2,label='t = {}s'.format(fmt(time)))

ax.set_xlabel('$p_f/p_0$')
ax.set_ylabel('z/L');
### ---------------------------------------------------------

# Numerical solution

path = 'gold/TerzaghiImportDataLayers_csv_pp_'
files = ['0001', '0054', '0043','0058']

for i, f in enumerate(files):
    df = pd.read_csv(path + f + '.csv')
    plt.plot( df.porepressure/p_0[0], df.y/100, 'k.', markersize=5, label='RHEA\'s solution' if i == 0 else '')

ax.legend(fontsize=20,ncol=2, handleheight=1.5, labelspacing=.1)
ax.grid(which = 'minor')

plt.savefig('Terzaghis_problem_layers.png', bbox_inches='tight')
