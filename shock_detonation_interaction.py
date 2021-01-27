# Transmitted Shockwave and Detonation Wave w/ Rarefaction Region
# Sam Faulk
""" with Perfect Gas Assumptions """
import sys
import os
os.chdir('/Users/samfaulk/Documents/git-sdf33/python-code')
print('Current Working Directory:',os.getcwd())
print('Currect Python Directory:',sys.executable)
import cantera as ct
from sdtoolbox import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from sympy import symbols , solve

img = mpimg.imread('detonation-shock-interaction.png')
imgplot = plt.imshow(img)
plt.show()

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
gasR = ct.Solution(mechR)
gasR.TPX = 1500.0 , 101325.0 , compR
GasR.P = gasR.P #pressure
GasR.T = gasR.T #temperature
GasR.M = gasR.mean_molecular_weight #molecular weight
GasR.rho = gasR.density #density
GasR.cp = gasR.cp_mass #specific heat, constant pressure
GasR.cv = gasR.cv_mass #specific heat, constant volume
GasR.hf = gasR.enthalpy_mole #enthalpy of formation

# Reactor
reactr = ct.IdealGasReactor(contents = gasR, energy = 'on')
reactr.volume = 10.0

# Simulation
simulation = ct.ReactorNet([reactr])

# Run the simulation
    # Assumed steady state of detonation is reached immediately
simulation.advance_to_steady_state()

# Product Gas after Simulation has run
GasP.T = gasR.T
GasP.P = gasR.P
GasP.M = gasR.mean_molecular_weight
GasP.X = gasR.X
GasP.rho = gasR.density
GasP.cp = gasR.cp_mass
GasP.cv = gasR.cv_mass
GasP.hf = gasR.enthalpy_mole

# Gas Constants
R.r = 1000*8.314/GasR.M
R.p = 1000*8.314/GasP.M

# Specific Heat Ratios
gamma = 1.2 # assume same gamma value

# Speed of Sound
GasP.c = np.sqrt(gamma*R.p*GasP.T)
GasR.c = np.sqrt(gamma*R.r*GasR.T)

# Mach of Incident Detonation Wave
def detMach_incident():
    [cj_speed] = CJspeed(GasP.P, GasP.T, GasP.X, mechP, fullOutput=True)
    Mdi = cj_speed/GasP.c
    return Mdi

# Mach of Incident Shock Wave
Msi = 5.0 # this is the incident shock Mach

def hf():
    delta_h = symbols('delta_h')
    Mdi = detMach_incident()
    expr = np.sqrt(((gamma**2.0-1.0)*delta_h)/(2.0*gamma*R.p*GasP.T)+1.0) + np.sqrt(((gamma**2.0-1.0)*delta_h)/(2.0*gamma*R.p*GasP.T)) - Mdi
    hf = solve(expr,delta_h)
    return hf

def postWave_incident():
    # (1): behind incident detonation wave
    # (2): between incident detonation and shock
    # (3): behind incident shock wave
    p1 = p2(gamma+1.0)/(gamma*Md**2.0+1.0) # pressure behind the incident detonation
    rho1 = rho2*(1.0+gamma*Md**2.0)/((gamma+1.0)*Md**2.0) # density behind the incident detonation
    p3 = GasR.P(gamma+1.0)/(2*gamma*(Ms**2.0-1.0)+(gamma+1.0)) # pressure behind the incident shock
    rho3 = GasR.P*((gamma-1.0)*Ms**2.0+2.0)/((gamma+1.0)*Ms**2.0) # density behind the incident shock
    Incident = np.array([p1,p3,rho1,rho3])
    return Incident

def Mach_transmit():
    Incident = postWave_incident()
    Mdi = detMach_incident()
    # Absolute Velocites in post wave regions, (1) & (3)
    w1 = Mdi*GasP.c
    w3 = Msi*GasR.c
    # Guess Mdt
    Mdt = np.arange(2.0,7.0,0.1)
    u_Mdi = Mdi*GasR.c + w3
    u_Msi = Msi*GasR.c + w1
    Mst = np.arange(0.1,25.0,0.1)
    # Properties from transmitted detonation equations
    p2_det = Incident[1]*(gamma+1.0)/(gamma*Mdt**2.0+1.0)
    rho2_det = Incident[3]*(1.0+gamma*Mdt**2.0)/((gamma+1.0)*Mdt**2.0)
    u2_det = u_Mdi*Incident[3]/rho2_det
    # Properties from transmitted shock equations
    p2_shock = Incident[0]*(gamma+1.0)/(2*gamma*(Mst**2.0-1.0)+(gamma+1.0))
    rho2_shock = Incident[2]*((gamma-1.0)*Mst**2.0+2.0)/((gamma+1.0)*Mst**2.0)
    u2_shock = u_Msi*Incident[2]/rho2_det

    print('p2_det: ',p2_det)
    print('p2_shock: ',p2_shock)
    print('\n')
    print('u2_det:',u2_det)
    print('u2_shock:',u2_shock)

output = Mach_transmit()







# Gas Constants
# R = R(1000*8.314/GasR.M , 1000*8.314/GasP.M)

# Specific Heat Ratios
# gamma = Gamma(GasR.cp/GasR.cv , GasP.cp/GasP.cv)

# Speed of Sound
# GasP.c = np.sqrt(gamma.p*R.p*GasP.T)
# GasR.c = np.sqrt(gamma.r*R.r*GasR.T)

# # Heat Released
# # Q = GasP.hf - GasR.hf
# # Q = 286000.0*2.0
# # gas_cj = PostShock_eq(cj_speed, GasP.P, GasP.T, compP, mechP)
# # gas_vn = PostShock_fr(cj_speed, GasP.P, GasP.T, compP, mechP)

# # Q = GasP.c*GasP.c*(Md-(1.0/Md))*(Md-(1.0/Md))/(2.0*((gamma**2.0) - 1.0))

# # k Value
# # k = 2.0*(((gamma.r)*(gamma.p-gamma.r)*(gamma.p+1.0)/((gamma.p**2.0)*(gamma.r-1.0))) + ((gamma.r**2.0)*(gamma.p**2.0 - 1.0)*(Q)/((gamma.p**2.0)*(GasR.c**2.0))))
# k = 2.0*(((gamma)*(gamma-gamma)*(gamma+1.0)/((gamma**2.0)*(gamma-1.0))) + ((gamma**2.0)*(gamma**2.0 - 1.0)*(Q)/((gamma**2.0)*(GasR.c**2.0))))

# # Detonation Strength
# # eta = Eta((gamma.r/gamma.p) + (k/2.0) - np.sqrt(k*((gamma.r/gamma.p)+(k/4.0))))
# eta.CJ = (gamma/gamma) + (k/2.0) - np.sqrt(k*((gamma/gamma)+(k/4.0)))

# def detMach():
#     # Products
#     p2 = GasP.P
#     rho2 = GasP.rho
#     # Reactants
#     p1 = GasR.P
#     rho1 = GasR.rho
#     hf = -241830.0*2.0 #https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI&Mask=1
#     Md = np.sqrt((p2*rho2)/(p1*rho1)*(2.0/(gamma.p-1.0) - 1.0) - (2.0*hf)/(gamma.p*p1*rho1) - (2.0)/(gamma.p-1.0))
#     return Md

# def detMach_incident():
#     # Products
#     p2 = GasP.P
#     rho2 = GasP.rho
#     # Reactants
#     p1 = GasR.P
#     rho1 = GasR.rho
#     hf = hf()
#     Md = np.sqrt((p2*rho2)/(p1*rho1)*(2.0/(gamma-1.0) - 1.0) - (2.0*hf)/(gamma*p1*rho1) - (2.0)/(gamma-1.0))
#     return Md

# def shockMach():
#     Mst = []
#     Ms = np.arange(5.0,24.0,1.0)
#     Md = detMach() #gives Md
#     print('Detonation Mach Number: %.3f'%Md)
#     index = np.arange(len(Ms))
#     # Products
#     c1 = GasP.c
#     p1 = GasP.P
#     # Reactants
#     c2 = GasR.c
#     p2 = GasR.P
#     expr1 = []
#     expr2 = []
#     for i in index:
#         u1 = GasP.c*Ms[i]
#         u2 = GasR.c*Md
#         expr1.append(u1 - c1*((2.0*gamma.p*Ms[i]**2.0-(gamma.p-1.0))/(gamma.p+1.0)-1.0)*np.sqrt(2.0/(gamma.p*(gamma.p-1.0)*(1.0+(2.0*gamma.p*Ms[i]**2.0-(gamma.p-1.0)/(gamma.p-1.0))))))
#         expr2.append(u2 - ((2.0*c2*gamma.p*(gamma.r+eta.CJ))/(np.sqrt(eta.CJ)*gamma.r*(gamma.p**2.0-1.0)))*(1.0 - ((p1)/(p2)*(eta.CJ*(2.0*gamma.p*Ms[i]**2.0-(gamma.p-1.0)))/(eta.CJ+gamma.r))**((gamma.p-1.0)/(2.0*gamma.p)) + ((gamma.p-1.0)*(gamma.p*eta.CJ-gamma.r))/((2.0*gamma.p)*(gamma.r+eta.CJ))))
#     exprs1 = [round(num1, 3) for num1 in expr1]
#     exprs2 = [round(num2, 3) for num2 in expr2]
#     print('EXPRESSION 1:\n',exprs1)
#     print('EXPRESSION 2:\n',exprs2)
#     for i in index:
#         for j in index:
#             if expr1[i] == expr2[j]:
#                 Mst.append(Ms)
#                 print('yes')
#     return Mst
# Mst = shockMach()
# print('Shock Mach Number: ',Mst)
  
# def shockMach():
#     Ms = 10.0
#     Md = 27.0
#     Mst = []
#     c1 = GasP.c
#     p1 = GasP.P
#     u1 = GasP.c*Ms

#     c2 = GasR.c
#     p2 = GasR.P
#     u2 = GasR.c*Md

#     expr1 = []
#     expr2 = []
#     expr1.append(u1 - c1*((2.0*gamma.p*Ms[i]**2.0-(gamma.p-1.0))/(gamma.p+1.0)-1.0)*np.sqrt(2.0/(gamma.p*(gamma.p-1.0)*(1.0+(2.0*gamma.p*Ms[i]**2.0-(gamma.p-1.0)/(gamma.p-1.0))))))
#     expr2.append(u2 - ((2.0*c2*gamma.p*(gamma.r+eta.CJ))/(np.sqrt(eta.CJ)*gamma.r*(gamma.p**2.0-1.0)))*(1.0 - ((p1)/(p2)*(eta.CJ*(2.0*gamma.p*Ms[i]**2.0-(gamma.p-1.0)))/(eta.CJ+gamma.r))**((gamma.p-1.0)/(2.0*gamma.p)) + ((gamma.p-1.0)*(gamma.p*eta.CJ-gamma.r))/((2.0*gamma.p)*(gamma.r+eta.CJ))))
#     print('expression1',expr1)
#     print('expression2',expr2)
#     for i in index:
#         for j in index:
#             if expr1[i] == expr2[j]:
#                 Mst.append(Ms)
#     return Mst
# shockMach()
# print(Mst)

