#!/usr/bin/env python3

# Solution to Terzaghi consolidation as presented in
# Section 2.2 of the online manuscript: Arnold Verruijt "Theory and Problems of Poroelasticity" Delft University of Technology 2013.  But note that the "sigma" in that paper is the negative of the stress in PorousFlow.

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Base properties
soil_height = 100
soil_bulk_modulus = 8.4E7
soil_shear_modulus = 6.25E7
fluid_bulk_modulus = 2.2E9
fluid_mobility = 1E-11 / 1E-3
soil_initial_porosity = 0.2
biot_coefficient = 0.9
normal_stress_on_top = 1000


# Derived properties
soil_lame_mu = soil_shear_modulus
soil_lame_lambda = soil_bulk_modulus - 2.0 * soil_shear_modulus / 3.0
soil_confined_compressibility = 1.0 / (soil_bulk_modulus + 4.0 * soil_lame_mu / 3.0)
soil_bulk_compliance = 1.0 / soil_bulk_modulus
fluid_bulk_compliance = 1.0 / fluid_bulk_modulus
soil_initial_storativity = soil_initial_porosity / fluid_bulk_modulus + (biot_coefficient - soil_initial_porosity) * (1.0 - biot_coefficient) / soil_bulk_modulus
consolidation_coefficient = fluid_mobility / (soil_initial_storativity + soil_confined_compressibility * biot_coefficient**2)
initial_porepressure = biot_coefficient * soil_confined_compressibility * normal_stress_on_top / (soil_initial_storativity + soil_confined_compressibility * biot_coefficient**2)
final_displacement = normal_stress_on_top * soil_height * soil_confined_compressibility
initial_displacement = normal_stress_on_top * soil_height * soil_confined_compressibility * soil_initial_storativity / (soil_initial_storativity + soil_confined_compressibility * biot_coefficient**2)

def expectedU(t):
    ctoverh2 = t * consolidation_coefficient / soil_height**2

    def coeff(n, ctoverh2):
        return np.power(1.0 / (2.0 * n - 1.0), 2) * np.exp(-np.power(2.0 * n - 1.0, 2) * np.pi ** 2 * ctoverh2 / 4.0)

    result = 0.0 * t
    for i in range(1, 11):
        result += coeff(i, ctoverh2)
    return 1.0 - (8.0 / np.pi**2) * result

def expectedP(t, z):
    ctoverh2 = t * consolidation_coefficient / soil_height**2
    zoverh = z / soil_height

    def coeff(n, zoverh):
        return np.power(-1.0, n - 1) / (2.0 * n - 1.0) * np.cos((2.0 * n - 1.0) * np.pi * zoverh / 2.0) * np.exp(-np.power(2 * n - 1, 2) * np.pi**2 * ctoverh2 / 4.0)

    result = 0.0 * z
    for i in range(1, 21):
        result += coeff(i, zoverh)
    approx = round(ctoverh2, -int(np.floor(np.sign(ctoverh2) * np.log10(abs(ctoverh2)))) + 1)
    return (approx, zoverh, (4.0 / np.pi) * result)

def get_moose_results(fn):
    f = open(fn, 'r')
    data = [list(map(float, line.strip().split(","))) for line in f.readlines()[1:] if line.strip()]
    f.close()
    return data

moose = get_moose_results("gold/TerzaghiImportData_full.csv")


plt.figure(1)
zpoints = np.arange(0, 100.1, 0.1)
zoverh = np.arange(0, 1.1, 0.1)
colors = ['k', 'k', 'r', 'b', 'g']

for i in range(1, len(moose)):
    t = moose[i][0]
    ex = expectedP(t, zpoints)
    plt.plot(ex[2], ex[1], colors[i] + '-', linewidth = 2.0, label = 'expected, t=' + str(t))
    plt.plot(np.array(moose[i][1:]) / initial_porepressure, zoverh, colors[i] + 'o', label = 'MOOSE, t=' + str(t))
plt.legend()
plt.xlabel("P/P0")
plt.ylabel("z/h")
plt.title("Terzaghi's consolidation problem")
plt.savefig("terzaghi_p.png")
sys.exit(0)
