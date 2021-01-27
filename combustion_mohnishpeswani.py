##Mohnish Peswani
##Constant volume process data is saved to csv file
##Script uses the curve fit result for Kef to return Kr1 and Kr2
##Output is saved to file
##Mixture, Temperature range and mechanism can be changed if necessary

import multiprocessing
import numpy as np
import cantera as ct
import itertools
from time import time
import pandas as pd
import csv
import matplotlib.pyplot as plt
from scipy.optimize import nnls
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import fsolve
import numpy as np
from scipy import interpolate, optimize
import sys
from scipy.signal import find_peaks

np.set_printoptions(threshold=sys.maxsize)

gases = {}


##Curve fit values changed here!
def kef_fit(T,rho):
    A= 1.0031540438356196e+28#1.78*10**27#6.667015647916756*10**26     #
    m= 0.9914410693955888#1.0#0.9936944229782382                      ##
    n=-4.180849809909199#-4.0 #-3.8258276679334515                       ##
    Ea= -229.22054296150841  #-224.3#-225.82556252393505
    kef = A*(rho**m)*(T**n)*np.exp(Ea*298.0/T)

    return kef

def init_process(mech):
    """
    This function is called once for each process in the Pool. We use it to
    initialize any Cantera objects we need to use.
    """
    gases[mech] = ct.Solution(mech)
    gases[mech].transport_model = 'Multi'

def run_simulation(args):
	tfinal = 100.0
	#tfinal = 1e-4
	tnow = 0.0
	told = 0.0
	mech,T,D,X = args
	Told=T
	gas = gases[mech]
	gas.TDX = T,D,X
	combustor = ct.IdealGasReactor(gas)
	sim=ct.ReactorNet([combustor])
	#sim.rtol = 1E-5
	#sim.atol = 1E-12
	gasR=ct.Solution('gri30.cti')
	gasP1=ct.Solution('gri30.cti')
	gasP2=ct.Solution('gri30.cti')

	csvfile='combustor_final.csv'
	outfile = open(csvfile,'w')
	csvwriter=csv.writer(outfile, lineterminator="\n")
	csvwriter.writerow(['#t (s)','T (K)','P (Pa)','rho (kg/s)']  + ['R','P1','P2'] + ['rhoP1','rhoP2','rhoR','WP1','WP2'] + gas.species_names)



	while tnow < tfinal:
		tnow = sim.step()
		#sim.set_max_time_step(1E-10)
		Told=combustor.T
		told=tnow
		H=combustor.thermo.enthalpy_mole
		W=combustor.thermo.mean_molecular_weight
		den=combustor.thermo.density

		gasR.TPX=combustor.T,combustor.thermo.P,'C2H2:2,O2:5'
		gasP1.TPX=combustor.T,combustor.thermo.P,'CO2:4,H2O:2'
		gasP2.TPX=combustor.T,combustor.thermo.P,'CO:4,H:4,O:6'
		WR=gasR.mean_molecular_weight
		print(WR)
	    #WR = gasR.mean_molecular_weight
		WP1 = gasP1.mean_molecular_weight
		WP2 = gasP2.mean_molecular_weight
		HR = gasR.enthalpy_mole
		HP1 = gasP1.enthalpy_mole
		HP2 = gasP2.enthalpy_mole

		matrix1= np.matrix([[WR,WP1,WP2],[HR,HP1,HP2],[1.0,1.0,1.0]])
		SolVec = ([W,H,1.0])

		moleFracs = matrix1.I.dot(SolVec)

		XR = (float(moleFracs[0,0]))
		XP1 = float(moleFracs[0,1])
		XP2 = float(moleFracs[0,2])

        #mass fractions:
		mass_frac_R = XR*(WR/W)
		mass_frac_P1=  XP1*(WP1/W)
		mass_frac_P2 = XP2*(WP2/W)

		rhoR = mass_frac_R*den
		rhoP1 = mass_frac_P1*den
		rhoP2 = mass_frac_P2*den


		csvwriter.writerow([tnow,  combustor.T, combustor.thermo.P, combustor.thermo.density]  + [XR,XP1,XP2] + [rhoP1,rhoP2,rhoR,WP1,WP2,combustor.thermo.mean_molecular_weight] + list(combustor.thermo.X) )
        outfile.close()
        ts,Temp,Pres,den,mFracR,mFracP1,mFracP2,rho1,rho2,rhoR,W = np.loadtxt('combustor_final.csv', delimiter=',', usecols=(0,1,2,3,4,5,6,7,8,9,12), unpack=True)

        result = ts,Temp,Pres,den,mFracR,mFracP1,mFracP2,rho1,rho2,rhoR,W

	return result

def serial(mech, predicate, nTemps):
    D0=100.0       #Tested in serial with slower output.
    D1=10.0
    D2=1.0
    D3=0.1
    D4=0.01
    D5=0.001

    T_initial= 1000.0
    T_final =  9000.0

    X = 'C2H2:2, O2:5'
    init_process(mech)


    y0 = list(map(predicate,
                 zip(itertools.repeat(mech),
                     np.linspace(T_initial, T_final, nTemps),
                     itertools.repeat(D0),
                     itertools.repeat(X))))

    y1 = list(map(predicate,
                 zip(itertools.repeat(mech),
                     np.linspace(T_initial, T_final, nTemps),
                     itertools.repeat(D1),
                     itertools.repeat(X))))


    y2 = list(map(predicate,
                 zip(itertools.repeat(mech),
                     np.linspace(T_initial, T_final, nTemps),
                     itertools.repeat(D2),
                     itertools.repeat(X))))

    y3 = list(map(predicate,
                 zip(itertools.repeat(mech),
                     np.linspace(T_initial, T_final, nTemps),
                     itertools.repeat(D3),
                     itertools.repeat(X))))


    y4 = list(map(predicate,
                 zip(itertools.repeat(mech),
                     np.linspace(T_initial, T_final,nTemps),
                     itertools.repeat(D4),
                     itertools.repeat(X))))



    y5 = list(map(predicate,
                 zip(itertools.repeat(mech),
                     np.linspace(T_initial, T_final, nTemps),
                     itertools.repeat(D5),
                     itertools.repeat(X))))
    return y0,y1,y2,y3,y4,y5

if __name__ == '__main__':
    nPoints = 15
    nProcs = 5

    print('Run in series')
    y0,y1,y2,y3,y4,y5 = serial('gri30.cti',run_simulation, nPoints)

ts_0,Temp_0,Pres_0,den_0,mFracR_0,mFracP1_0,mFracP2_0,rho1_0,rho2_0,rhoR_0,W_0= map(list, zip(*y0))
ts_1,Temp_1,Pres_1,den_1,mFracR_1,mFracP1_1,mFracP2_1,rho1_1,rho2_1,rhoR_1,W_1= map(list, zip(*y1))
ts_2,Temp_2,Pres_2,den_2,mFracR_2,mFracP1_2,mFracP2_2,rho1_2,rho2_2,rhoR_2,W_2= map(list, zip(*y2))
ts_3,Temp_3,Pres_3,den_3,mFracR_3,mFracP1_3,mFracP2_3,rho1_3,rho2_3,rhoR_3,W_3= map(list, zip(*y3))
ts_4,Temp_4,Pres_4,den_4,mFracR_4,mFracP1_4,mFracP2_4,rho1_4,rho2_4,rhoR_4,W_4= map(list, zip(*y4))
ts_5,Temp_5,Pres_5,den_5,mFracR_5,mFracP1_5,mFracP2_5,rho1_5,rho2_5,rhoR_5,W_5= map(list, zip(*y5))

f= open("combustor_data_konnov_kr1_5","w")
np.savez(f,ts_0=ts_0,Temp_0=Temp_0,Pres_0=Pres_0,den_0=den_0,mFracR_0=mFracR_0,mFracP1_0=mFracP1_0,mFracP2_0=mFracP2_0,rho1_0=rho1_0,rho2_0=rho2_0,rhoR_0=rhoR_0,W_0=W_0,
          ts_1=ts_1,Temp_1=Temp_1,Pres_1=Pres_1,den_1=den_1,mFracR_1=mFracR_1,mFracP1_1=mFracP1_1,mFracP2_1=mFracP2_1,rho1_1=rho1_1,rho2_1=rho2_1,rhoR_1=rhoR_1,W_1=W_1,
          ts_2=ts_2,Temp_2=Temp_2,Pres_2=Pres_2,den_2=den_2,mFracR_2=mFracR_2,mFracP1_2=mFracP1_2,mFracP2_2=mFracP2_2,rho1_2=rho1_2,rho2_2=rho2_2,rhoR_2=rhoR_2,W_2=W_2,
          ts_3=ts_3,Temp_3=Temp_3,Pres_3=Pres_3,den_3=den_3,mFracR_3=mFracR_3,mFracP1_3=mFracP1_3,mFracP2_3=mFracP2_3,rho1_3=rho1_3,rho2_3=rho2_3,rhoR_3=rhoR_3,W_3=W_3,
          ts_4=ts_4,Temp_4=Temp_4,Pres_4=Pres_4,den_4=den_4,mFracR_4=mFracR_4,mFracP1_4=mFracP1_4,mFracP2_4=mFracP2_4,rho1_4=rho1_4,rho2_4=rho2_4,rhoR_4=rhoR_4,W_4=W_4,
          ts_5=ts_5,Temp_5=Temp_5,Pres_5=Pres_5,den_5=den_5,mFracR_5=mFracR_5,mFracP1_5=mFracP1_5,mFracP2_5=mFracP2_5,rho1_5=rho1_5,rho2_5=rho2_5,rhoR_5=rhoR_5,W_5=W_5)