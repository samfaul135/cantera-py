# Viscosities & Prandtl Number
# function with density input
# Sam Faulk
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import time
import math
# Function Density
def Density():
    # Set Limits for Density
    dens_initial = 0.001
    dens_now = dens_initial
    dens_final = 100.0
    dens_delta = 10.0
    dens_mag = []
    while dens_now <= dens_final:
        dens_mag.append(dens_now)
        dens_now *= dens_delta
    # Create class Visc
    class Visc:
        pass
    visc = Visc()
    # Set Mechanism
    mech = 'gri30.cti'
    mixing = 'gri30_mix'
    # Set Composition
    comp = 'O2:2.0 , CH4:1.0'
    # Set Limits of Temperature
    temp_initial = 300.0
    temp_final = 3000.0
    temp_delta = 500.0
    temp = np.arange(temp_initial,temp_final,temp_delta)
    # Create Indexing Arrays
    index = np.arange(len(dens_mag))
    element = np.arange(len(temp))
    # Create Viscosity Array
    visc.d = np.zeros((len(dens_mag),len(temp)))
    visc.k = np.zeros((len(dens_mag),len(temp)))
    # Create Prandtl Number Array
    Pr = np.zeros((len(dens_mag),len(temp)))
    for i in index: # for density
        print(i)
        for e in element: # for temperature
            print(e)
            # Create gasA
            gasA = ct.Solution(mech,mixing)
            gasA.TDX = temp[e] , dens_mag[i] , comp
            # Dynamic Viscosity
            visc.d[i,e] = gasA.viscosity
            # Kinematic Viscosity
            visc.k[i,e] = gasA.viscosity/gasA.density
            # Prandtl Number
            Pr[i,e] = gasA.cp*gasA.viscosity/gasA.thermal_conductivity
    # Create temp array for 
    mag = np.log10(dens_final / dens_initial) # define increment of change to ensure that temp and dens are the same size
    temp_delta = (temp_final - temp_initial) / mag # define NEW delta increment so that "dens" and "temp" lists are same size for plotting
    temp_array = np.arange(temp_initial , temp_final , temp_delta)
    # Create class States()
    class States():
        pass
    # Setting values of the class
    states = States()
    states.viscd = visc.d
    states.visck = visc.k
    states.Pr = Pr
    states.index = index
    states.element = element
    states.temp_array = temp_array
    states.temp_final = temp_final
    return states
# Call function
states = Density()
# Append final temperature onto end of temp array
states.temp_array = np.append(states.temp_array , states.temp_final)
# Plot Dynamic Viscosity
plt.figure(1)
for i in states.index: # for density
    plt.plot(states.temp_array, states.viscd[i,:])
plt.xlabel('$Initial$ $Temperature$ (K)')
plt.ylabel('$Dynamic$ $Viscosity$ (Pa-s)')
plt.legend(('0.001 [kg/m^3]' , '0.01' , '0.1' , '1.0' , '10.0' , '100.0') , loc = 'best')
plt.grid()
plt.tight_layout()
plt.show()
# Plot Kinematic Viscosity
plt.figure(2)
for i in states.index: # for density
    plt.plot(states.temp_array, states.visck[i,:])
plt.xlabel('$Initial$ $Temperature$ (K)')
plt.ylabel('$Kinematic$ $Viscosity$ (m^2/s)')
plt.legend(('0.001 [kg/m^3]' , '0.01' , '0.1' , '1.0' , '10.0' , '100.0') , loc = 'best')
plt.grid()
plt.tight_layout()
plt.show()
# Plot Prandtl Number
plt.figure(3)
for i in states.index: # for density
    plt.plot(states.temp_array, states.Pr[i,:])
plt.xlabel('$Initial$ $Temperature$ (K)')
plt.ylabel('$Prandtl$ $Number$')
plt.legend(('0.001 [kg/m^3]' , '0.01' , '0.1' , '1.0' , '10.0' , '100.0') , loc = 'best')
plt.grid()
plt.tight_layout()
plt.show()