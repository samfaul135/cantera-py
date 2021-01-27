# Viscosities & Prandtl Number
# Sam Faulk
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import time
import math
dVisc = []
kVisc = []
Pr = []
temp = []
dens = []
# Set limits for density
dens_initial = 0.001
dens_now = dens_initial
dens_final = 100.0
dens_delta = 10.0
# Create increment of change to ensure that temp and dens are the same size
mag = np.log10(dens_final / dens_initial)
while dens_now <= dens_final:
    # Keeps track of progress, as increments of temperture
    print('{}kg/m^3 of {}kg/m^3'.format(dens_now,dens_final))
    # Set limits for temperature
    temp_initial = 300.0
    temp_now = temp_initial
    temp_final = 3000.0
    temp_delta = (temp_final - temp_initial) / mag # define delta increment that "dens" and "temp" lists are same size
    while temp_now <= temp_final:
        print('    {}K of {}K'.format(temp_now,temp_final))
        # Define empty lists for within the "while loop"
        dVisc_within = []
        kVisc_within = []
        Pr_within = []
        # Create gasA
        gasA = ct.Solution('gri30.cti','gri30_mix')
        gasA.TDX = temp_now , dens_now , 'O2:2.0 , CH4:1.0'
        # Dynamic Viscosity
        dVisc_within.append(gasA.viscosity)
        dVisc = np.append(dVisc , dVisc_within , axis = 0)
        # Kinematic Viscosity
        kVisc_within.append(gasA.viscosity/gasA.density)
        kVisc = np.append(kVisc , kVisc_within , axis = 0)
        # Prandtl Number
        Pr_within.append(gasA.cp*gasA.viscosity/gasA.thermal_conductivity)
        Pr = np.append(Pr , Pr_within , axis = 0)
        # Temperature increments
        temp_now += temp_delta
    # Density list
    dens.append(dens_now)
    # Density increments
    dens_now *= dens_delta
# Temperature list
temp = np.arange(temp_initial , temp_final , temp_delta)
temp = np.append(temp , temp_final)
# Split Viscosity lists
mag = int(mag)
dVisc = np.split(dVisc,mag+1)
kVisc = np.split(kVisc,mag+1)
Pr = np.split(Pr,mag+1)
# Plot values
rnge = np.arange(mag+1)
# Plot Dynamic Viscosity
plt.figure(1)
for i in rnge:
    plt.plot(temp, dVisc[i])
plt.xlabel('$Initial$ $Temperature$ (K)')
plt.ylabel('$Dynamic$ $Viscosity$ (Pa-s)')
plt.legend(('0.001 [kg/m^3]' , '0.01' , '0.1' , '1.0' , '10.0' , '100.0') , loc = 'best')
plt.grid()
plt.tight_layout()
plt.show()
# Plot Kinematic Viscosity
plt.figure(2)
for i in rnge:
    plt.plot(temp, kVisc[i])
plt.xlabel('$Initial$ $Temperature$ (K)')
plt.ylabel('$Kinematic$ $Viscosity$ (m^2/s)')
plt.legend(('0.001 [kg/m^3]' , '0.01' , '0.1' , '1.0' , '10.0' , '100.0') , loc = 'best')
plt.grid()
plt.tight_layout()
plt.show()
# Plot Prandtl Number
plt.figure(3)
for i in rnge:
    plt.plot(temp, Pr[i])
plt.xlabel('$Initial$ $Temperature$ (K)')
plt.ylabel('$Prandtl$ $Number$')
plt.legend(('0.001 [kg/m^3]' , '0.01' , '0.1' , '1.0' , '10.0' , '100.0') , loc = 'best')
plt.grid()
plt.tight_layout()
plt.show()