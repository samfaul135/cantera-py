# Cantera Reactor Basics
# Sam Faulk
# Modeling an Ideal Gas H2 Combustion
import cantera as ct
import matplotlib.pyplot as plt
import time

# Time Start
#tic = time.clock()

# Define inital conditions
T_0 = 1500.0 #K
P_0 = 101325.0 #Pa

# Import a gas mechanism
mech = 'gri30.cti'
gas1 = ct.Solution(mech)

# Define state properties of the gas
comp = 'H2:2.0 , O2:1.0, N2:3.76'
gas1.TPX = T_0 , P_0 , comp

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
t_final = 1e-4

while t_now < t_final:
    temp.append(reactr.T)
    press.append(reactr.thermo.P)
    states.append(reactr.thermo.state)
    time.append(t_now)
    t_now = simulation.step()

p_max = max(press)/101325.0
print('Constant Volume Ignition of ' +comp+ ' using the ' +mech+ ' mechanism')
print('Maximum Pressure : %.4f atm'%p_max)

plt.subplot(2, 2, 1)
plt.plot(time, temp)
plt.xlabel('$Time$ (sec)')
plt.ylabel('$Temperature$ (K)')
plt.grid()

plt.subplot(2, 2, 2)
plt.plot(time, press)
plt.xlabel('$Time$ (sec)')
plt.ylabel('$Pressure$ (Pa)')
plt.grid()

plt.subplot(2, 2, 3)
plt.plot(time, states.X[:,gas1.species_index('H2')])
plt.xlabel('$Time$ (sec)')
plt.ylabel('$H2$ Mole Fraction')
plt.grid()

plt.subplot(2, 2, 4)
plt.plot(time, states.X[:,gas1.species_index('H2O')])
plt.xlabel('$Time$ (sec)')
plt.ylabel('$H2O$ Mole Fraction')
plt.grid()

plt.tight_layout()
plt.show()

simulation.advance_to_steady_state()

print(reactr.thermo.P)