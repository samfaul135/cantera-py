# Transmitted Shockwave and Detonation Wave
# Sam Faulk
""" with Perfect Gas Assumptions """
#import cantera as ct
#from sdtoolbox.postshock import CJspeed
#from sdtoolbox.postshock import PostShock_eq
#from sdtoolbox.thermo import soundspeed_eq, soundspeed_fr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


# img = mpimg.imread('detonation-shock-interaction.png')
# imgplot = plt.imshow(img)
# plt.show()

# Set Mechanism
mechR = 'gri30.cti'
mechP = 'gri30.cti'

# Set Composition
compR = 'H2:2.0, O2:1.0, N2:3.76'


class Gamma:
    pass


class Eta:
    pass


class GasR:
    pass


class GasP:
    pass


class R:
    pass


# Defining Classes
eta = Eta()
GasR = GasR()
GasP = GasP()
R = R()

# Reactant Gas
#gasR = ct.Solution(mechR)
#gasR.TPX = 1500.0, 101325.0, compR
# GasR.P = gasR.P  # pressure
GasR.P = 101325.0
# GasR.T = gasR.T  # temperature
GasR.T = 1500.0
# GasR.M = gasR.mean_molecular_weight  # molecular weight
# GasR.rho = gasR.density  # density
GasR.rho = 0.16989275549920685
# GasR.cp = gasR.cp_mass  # specific heat, constant pressure
# GasR.cv = gasR.cv_mass  # specific heat, constant volume
# GasR.hf = gasR.enthalpy_mole  # enthalpy of formation

# Reactor
#reactr = ct.IdealGasReactor(contents=gasR, energy='on')
#reactr.volume = 10.0

# Simulation
#simulation = ct.ReactorNet([reactr])

# Run the simulation
# Assumed steady state of detonation is reached immediately
# simulation.advance_to_steady_state()

# Product Gas after Simulation has run
#GasP.T = gasR.T
#GasP.P = gasR.P
GasP.P = 183658.94369651366
#GasP.M = gasR.mean_molecular_weight
#GasP.X = gasR.X
#GasP.rho = gasR.density
GasP.rho = 0.16989275549920685
#GasP.cp = gasR.cp_mass
#GasP.cv = gasR.cv_mass
#GasP.hf = gasR.enthalpy_mole

# Gas Constants
R.r = 1000*8.314/GasR.M
R.p = 1000*8.314/GasP.M

# Specific Heat Ratios
gamma = 1.2  # assume same gamma value

# Speed of Sound
#GasP.c = np.sqrt(gamma*R.p*GasP.T)
#GasR.c = np.sqrt(gamma*R.r*GasR.T)
GasP.c = 1138.930347294427
GasR.c = 845.9592240707273


# def detMach_incident():  # Mach of Incident Detonation Wave
#   [cj_speed, R2, plot_data] = CJspeed(
#      GasR.P, GasR.T, compR, mechR, fullOutput=True)
# Mdi = cj_speed/GasR.c
# print(cj_speed)
# print(GasR.c)
# print(Mdi)
#  return Mdi


#Mdi = detMach_incident()
Mdi = 2.1463872338905077
# def hf():
#     #delta_h = symbols('delta_h')
#     Mdi = detMach_incident()
#     # expr = np.sqrt(((gamma**2.0-1.0)*delta_h)/(2.0*gamma*R.p*GasP.T)+1.0) + \
#     #    np.sqrt(((gamma**2.0-1.0)*delta_h)/(2.0*gamma*R.p*GasP.T)) - Mdi

#     #hf = solve(expr, delta_h)
#     hf = (gamma/(gamma-1.0))*GasP.P*GasP.rho*(1.0+(gamma-1.0)/2.0) - \
#         (gamma/(gamma-1.0))*GasR.P*GasR.rho*(1.0+Mdi**(2.0)*((gamma-1.0)/2.0))
#     return hf


# def postWave_incident():

#     # (1): behind incident detonation wave
#     # (2): between incident detonation and shock
#     # (3): behind incident shock wave
#     # pressure behind the incident detonation

#     # Mach of Incident Shock Wave
#     Msi = 5.0  # this is the incident shock Mach
#     Mdi = detMach_incident()
#     p2 = GasR.P
#     rho2 = GasR.rho
#     p1 = p2*(gamma*Mdi**2.0+1.0)/(gamma+1.0)
#     # density behind the incident detonation
#     rho1 = rho2*((gamma+1.0)*Mdi**2.0)/(1.0+gamma*Mdi**2.0)
#     # pressure behind the incident shock
#     p3 = p2*(2*gamma*(Msi**2.0-1.0)+(gamma+1.0))/(gamma+1.0)
#     rho3 = rho2*((gamma+1.0)*Msi ** 2.0)/((gamma-1.0)*Msi **
#                                           2.0+2.0)  # density behind the incident shock
#     return p1, p3, rho1, rho3


# def Mach_transmit():
#     p1, p3, rho1, rho3 = postWave_incident()
#     Mdi = detMach_incident()

#     # Mach of Incident Shock Wave
#     Msi = 5.0  # this is the incident shock Mach
#     # Absolute Velocites in post wave regions, (1) & (3)
#     w1 = Mdi*GasP.c
#     w3 = Msi*GasR.c
#     # Guess Mdt and Mst
#     #Mdt = np.arange(2.0, 7.0, 0.1)
#     Mdt = np.arange(0.1, 25.0, 0.1)
#     Mst = np.arange(0.1, 25.0, 0.1)
#     u_Mdt = Mdt*GasP.c + w3  # velocities of wave with respect to lab frame
#     u_Mst = Mst*GasP.c + w1  # velocities of wave with respect to lab frame
#     # Properties from transmitted detonation equations
#     p2_det = p3*(gamma+1.0)/(gamma*Mdt**2.0+1.0)
#     rho2_det = rho3*(1.0+gamma*Mdt**2.0)/((gamma+1.0)*Mdt**2.0)
#     u2_det = u_Mdt*rho3/rho2_det
#     # Properties from transmitted shock equations
#     p2_shock = p1*(gamma+1.0)/(2*gamma*(Mst**2.0-1.0)+(gamma+1.0))
#     rho2_shock = rho1*((gamma-1.0)*Mst**2.0+2.0)/((gamma+1.0)*Mst**2.0)
#     u2_shock = u_Mst*rho3/rho2_shock

#     # print('p2_det: ', p2_det)
#     # print('p2_shock: ', p2_shock)
#     # print('\n')
#     # print('u2_det:', u2_det)
#     # print('u2_shock:', u2_shock)

#     tol = 1.0
#     press_det = []
#     press_shock = []
#     vel_det = []
#     vel_shock = []

#     for det in p2_det:
#         for shock in p2_shock:
#             if np.abs(det-shock) <= tol:
#                 #print('pressure match')
#                 press_det.append(Mdt[np.where(p2_det == det)])
#                 press_shock.append(Mst[np.where(p2_shock == shock)])

#     # detonation & shock Machs according to matching pressures
#     print('Mach Numbers calculated by matching pressures at contact surface')
#     print('     Mdt:%.3f' % np.mean(press_det))
#     print('     Mst:%.3f' % np.mean(press_shock))

#     for det in u2_det:
#         for shock in u2_shock:
#             if np.abs(det-shock) <= tol:
#                 #print("velocity match")
#                 vel_det.append(Mdt[np.where(u2_det == det)])
#                 vel_shock.append(Mst[np.where(u2_shock == shock)])
#     # detonation & shock Machs according to matching velocities
#     print('Mach Numbers calculated by matching velocities at contact surface')
#     print('     Mdt:%.3f' % np.mean(vel_det))
#     print('     Mst:%.3f' % np.mean(vel_shock))


# output = Mach_transmit()
