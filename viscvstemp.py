# Viscosity vs Temperature 
# Using various theories of gases
import numpy as np
import scipy as sp
from scipy.constants import k
import matplotlib.pyplot as plt

T = np.arange(300.0,3000.0,1.0) #K

# Hard-sphere kinetic theory
CH4_moles = 1.0 #mole
CH4_molarmass = 16.04 #g/mole
m = CH4_molarmass*CH4_moles/1000 #kg
kB = k #Boltzmann Constant
d = 1.0 #m
mu = 1.016*5.0*(16.0*d**2.0)*np.sqrt(k*m*T/np.pi)

plt.figure(1)
plt.plot(T,mu,label='Hard-sphere kinetic theory')
plt.xlabel('Temperature (K)')
plt.ylabel('Dynamic Viscosity (Pa-s)')
plt.legend(loc='best')

# Sutherland model
S = 197.8 # http://zigherzog.net/stirling/simulations/DHT/ViscosityTemperatureSutherland.html
T_ref = 300.0 #K
mu_ref = 11.13e-6 #Pa-s; https://www.engineeringtoolbox.com/methane-dynamic-kinematic-viscosity-temperature-pressure-d_2068.html?vA=300&degree=K#
mu = 5.0*(16.0*d**2.0)*np.sqrt(kB*m*T/np.pi)*(1.0+S/T)**(-1.0)
plt.figure(2)
plt.plot(T,mu,label='Sutherland model')

mu = mu_ref*(T/T_ref)**(1.5)*((T_ref+S)/(T+S))
plt.plot(T,mu,label='Sutherland model // from reference')
plt.xlabel('Temperature (K)')
plt.ylabel('Dynamic Viscosity (Pa-s)')
plt.legend(loc='best')

# Power Law Force
s = 0.5
T_prime = 300.0 #K
mu_prime = 11.13e-6 #dynamic viscosity, Pa-s; https://www.engineeringtoolbox.com/methane-dynamic-kinematic-viscosity-temperature-pressure-d_2068.html?vA=300&degree=K#
mu = mu_prime*(T/T_prime)**s
plt.figure(3)
plt.plot(T,mu,label='Power Law Force')
plt.xlabel('Temperature (K)')
plt.ylabel('Dynamic Viscosity (Pa-s)')
plt.legend(loc='best')


# Lennard-Jones
ep = 154.0*kB # for methane; Bird, Stewart, & Lightfoot (2007), pp. 864–865
Tstar = T*kB/ep
omega = 1.16145*Tstar**(-0.14874) + 0.52487*np.exp(-0.773207*Tstar) + 2.16178*np.exp(-2.437877*Tstar) # Neufeld, Jansen, & Aziz (1972)
sigma = 3.780 # for methane; Bird, Stewart, & Lightfoot (2007), pp. 864–865
mu = 5.0/(16.0*np.sqrt(np.pi))*np.sqrt(m*kB*T)/(sigma**(2.0)*omega)
plt.figure(4)
plt.plot(T,mu,label='Lennard-Jones')
plt.xlabel('Temperature (K)')
plt.ylabel('Dynamic Viscosity (Pa-s)')
plt.legend(loc='best')
plt.show()
