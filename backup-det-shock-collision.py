# Transmitted Shockwave and Detonation Wave w/ Rarefaction Region
# Sam Faulk
""" with Perfect Gas Assumptions """
import cantera as ct
from sdtoolbox.postshock import CJspeed
from sdtoolbox.postshock import PostShock_eq
from sdtoolbox.thermo import soundspeed_eq, soundspeed_fr
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
gasR = ct.Solution(mechR)
gasR.TPX = 1500.0, 101325.0, compR
GasR.P = gasR.P  # pressure
GasR.T = gasR.T  # temperature
GasR.M = gasR.mean_molecular_weight  # molecular weight
GasR.rho = gasR.density  # density
GasR.cp = gasR.cp_mass  # specific heat, constant pressure
GasR.cv = gasR.cv_mass  # specific heat, constant volume
GasR.hf = gasR.enthalpy_mole  # enthalpy of formation

# Reactor
reactr = ct.IdealGasReactor(contents=gasR, energy='on')
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
gamma = 1.2  # assume same gamma value

# Speed of Sound
GasP.c = np.sqrt(gamma*R.p*GasP.T)
GasR.c = np.sqrt(gamma*R.r*GasR.T)


def detMach_incident():  # Mach of Incident Detonation Wave
    [cj_speed, R2, plot_data] = CJspeed(
        GasR.P, GasR.T, compR, mechR, fullOutput=True)
    Mdi = cj_speed/GasR.c
    # print(Mdi)
    return Mdi


def hf():
    #delta_h = symbols('delta_h')
    Mdi = detMach_incident()
    # expr = np.sqrt(((gamma**2.0-1.0)*delta_h)/(2.0*gamma*R.p*GasP.T)+1.0) + \
    #    np.sqrt(((gamma**2.0-1.0)*delta_h)/(2.0*gamma*R.p*GasP.T)) - Mdi

    #hf = solve(expr, delta_h)
    hf = (gamma/(gamma-1.0))*GasP.P*GasP.rho*(1.0+(gamma-1.0)/2.0) - \
        (gamma/(gamma-1.0))*GasR.P*GasR.rho*(1.0+Mdi**(2.0)*((gamma-1.0)/2.0))
    return hf


def postWave_incident():

    # (1): behind incident detonation wave
    # (2): between incident detonation and shock
    # (3): behind incident shock wave
    # pressure behind the incident detonation

    # Mach of Incident Shock Wave
    Msi = 5.0  # this is the incident shock Mach
    Mdi = detMach_incident()
    p2 = GasR.P
    rho2 = GasR.rho
    p1 = p2*(gamma+1.0)/(gamma*Mdi**2.0+1.0)
    # density behind the incident detonation
    rho1 = rho2*(1.0+gamma*Mdi**2.0)/((gamma+1.0)*Mdi**2.0)
    # pressure behind the incident shock
    p3 = p2*(gamma+1.0)/(2*gamma*(Msi**2.0-1.0)+(gamma+1.0))
    rho3 = rho2*((gamma-1.0)*Msi**2.0+2.0)/((gamma+1.0)*Msi **
                                            2.0)  # density behind the incident shock
    return p1, p3, rho1, rho3


def Mach_transmit():
    p1, p3, rho1, rho3 = postWave_incident()
    Mdi = detMach_incident()

    # Mach of Incident Shock Wave
    Msi = 5.0  # this is the incident shock Mach
    # Absolute Velocites in post wave regions, (1) & (3)
    w1 = Mdi*GasP.c  # velocity of gases behind detonation wave
    w3 = Msi*GasR.c  # velocity of shocked gases
    # Guess Mdt and Mst
    #Mdt = np.arange(2.0, 7.0, 0.1)
    Mdt = np.arange(0.1, 25.0, 0.1)
    Mst = np.arange(0.1, 25.0, 0.1)
    u_Mdt = Mdt*GasP.c + w3  # velocities of wave with respect to lab frame
    u_Mst = Mst*GasP.c + w1  # velocities of wave with respect to lab frame
    # Properties from transmitted detonation equations
    p2_det = p3*(gamma+1.0)/(gamma*Mdt**2.0+1.0)
    rho2_det = rho3*(1.0+gamma*Mdt**2.0)/((gamma+1.0)*Mdt**2.0)
    u2_det = u_Mdt*rho3/rho2_det
    # Properties from transmitted shock equations
    p2_shock = p1*(gamma+1.0)/(2*gamma*(Mst**2.0-1.0)+(gamma+1.0))
    rho2_shock = rho1*((gamma-1.0)*Mst**2.0+2.0)/((gamma+1.0)*Mst**2.0)
    u2_shock = u_Mst*rho3/rho2_shock

    # print('p2_det: ', p2_det)
    # print('p2_shock: ', p2_shock)
    # print('\n')
    # print('u2_det:', u2_det)
    # print('u2_shock:', u2_shock)

    tol = 1.0
    press_det = []
    press_shock = []
    vel_det = []
    vel_shock = []

    for det in p2_det:
        for shock in p2_shock:
            if np.abs(det-shock) <= tol:
                #print('pressure match')
                press_det.append(Mdt[np.where(p2_det == det)])
                press_shock.append(Mst[np.where(p2_shock == shock)])

    # detonation & shock Machs according to matching pressures
    print('Mdt:%.3f' % np.mean(press_det))
    print('Mst:%.3f' % np.mean(press_shock))

    for det in u2_det:
        for shock in u2_shock:
            if np.abs(det-shock) <= tol:
                #print("velocity match")
                vel_det.append(Mdt[np.where(u2_det == det)])
                vel_shock.append(Mst[np.where(u2_shock == shock)])
    # detonation & shock Machs according to matching velocities
    print('Mdt:%.3f' % np.mean(vel_det))
    print('Mst:%.3f' % np.mean(vel_shock))


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
