## Cantera Reactor Basics
# Sam Faulk
# Modeling an Ideal Gas H2 Combustion
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# Define inital conditions
T_0 = 1500.0 #K
P_0 = 101325 #Pa

# Import a gas mechanism
gas1 = ct.Solution('gri30.cti')

# Define state properties of the gas
gas1.TPX = T_0 , P_0 , 'H2:2.0 , O2:1.0, N2:3.76'

# Create the reactor
reactr = ct.IdealGasReactor(contents = gas1, energy = 'on')
reactr.volume = 1.0

# Now create a reactor network
    # Note: this script only contains one reactor so only one is defined in ReactorNet([])
simulation = ct.ReactorNet([reactr])

# Run the simulation over decreasing time intervals until the reactor is extinguished
states = ct.SolutionArray(gas1) 
time = []
press = []
temp = []
t_initial = 0.0
t_now = 0.0
t_final = 1e-2

while t_now < t_final:
    temp.append(reactr.T)
    press.append(reactr.thermo.P)
    states.append(reactr.thermo.state)
    time.append(t_now)
    t_now = simulation.step()

p_max = max(press)/101325.0
print('Maximum Pressure of Constant Volume Ignition : %f atm' %p_max)

plt.subplot(2, 2, 1)
plt.plot(time, temp)
plt.xlabel('$Time$ (ms)')
plt.ylabel('$Temperature$ (K)')

plt.subplot(2, 2, 2)
plt.plot(time, press)
plt.xlabel('$Time$ (ms)')
plt.ylabel('$Pressure$ (Pa)')

plt.subplot(2, 2, 3)
plt.plot(time, states.X[:,gas1.species_index('H2')])
plt.xlabel('$Time$ (ms)')
plt.ylabel('$H2$ Mole Fraction')

plt.subplot(2, 2, 4)
plt.plot(time, states.X[:,gas1.species_index('H2O')])
plt.xlabel('$Time$ (ms)')
plt.ylabel('$H2O$ Mole Fraction')

plt.tight_layout()
plt.show()