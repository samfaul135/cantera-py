# Diffusion Coefficients, Viscosities, & Prandtl Number
# Sam Faulk
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import time
import math

tic = time.clock()

visc_d = []
visc_k = []
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
        visc_d_within = []
        visc_k_within = []

        # Create gasA
        gasA = ct.Solution('gri30.cti','gri30_mix')
        gasA.TPX = temp_now , 101325.0 , 'O2:0.21 , N2:0.78 , AR:0.01'
        #densA = gasA.density

        # Create gasB
        gasB = ct.Solution('gri30.cti','gri30_mix')
        gasB.TDX = temp_now , dens_now , 'CH4:1.0'

        # Create reserviorA, ReserviorB, and downstream
        #resA = ct.Reservoir(gasA)
        resB = ct.Reservoir(gasB)
        downstream = ct.Reservoir(gasA)

        # Redefine contents of gasB to put in reactor
        #gasB.TPX = temp_now , 101325.0 , 'O2:0.21 , N2:0.78 , AR:0.01'
        mix = ct.IdealGasReactor(gasA)

        # Create massflow controllers
        #mfc1 = ct.MassFlowController(resA , mix , mdot = densA*2.5/0.21)
        mfc2 = ct.MassFlowController(resB , mix , mdot = dens_now*1.0)

        # Connect 'mix' to 'downstream' with a valve
        outlet = ct.Valve(mix , downstream, K = 10.0)

        # Create reactor network 
        simul = ct.ReactorNet([mix])

        # Move simulation up to equilibrium, with constant enthalpy and pressure
        #gasB.equilibrate('HP',"auto")

        # Dynamic Viscosity
        visc_d_within.append(gasB.viscosity)
        visc_d = np.append(visc_d , visc_d_within , axis = 0)

        # Kinematic Viscosity
        visc_k_within.append(gasB.viscosity/dens_now)
        visc_k = np.append(visc_k , visc_k_within , axis = 0)

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
visc_d = np.split(visc_d,mag+1)
visc_k = np.split(visc_k,mag+1)

# Print time elapsed
toc = time.clock()
print('seconds:' +str(toc-tic))

# Plot values
rnge = np.arange(mag+1)

# Plot Dynamic Viscosity
plt.figure(1)
for i in rnge:
    plt.plot(temp, visc_d[i])
plt.xlabel('$Initial$ $Temperature$ (K)')
plt.ylabel('$Dynamic$ $Viscosity$ (Pa-s)')
plt.legend(('0.001 [kg/m^3]' , '0.01' , '0.1' , '1.0' , '10.0' , '100.0') , loc = 'best')
plt.grid()
plt.tight_layout()
plt.show()

# Plot Kinematic Viscosity
plt.figure(3)
for i in rnge:
    plt.plot(temp, visc_k[i])
plt.xlabel('$Initial$ $Temperature$ (K)')
plt.ylabel('$Kinematic$ $Viscosity$ (m^2/s)')
plt.legend(('0.001 [kg/m^3]' , '0.01' , '0.1' , '1.0' , '10.0' , '100.0') , loc = 'best')
plt.grid()
plt.tight_layout()
plt.show()