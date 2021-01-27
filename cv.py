"""
Constant volume ignition problems

December 2015
"""

import cantera as ct
import csv
import math
from optparse import OptionParser
parser = OptionParser()
OneAtm = 101300	# Pa
parser.add_option("--r", type="float", dest="den", default=1.0)
parser.add_option("--te", type="float", dest="te", default=1000.0)
parser.add_option("--x", type="string", dest="comp", default='CH4:1,O2:2')
parser.add_option("--mech", type="string", dest="mech", default='gri30_highT.cti')
(init_state, args) = parser.parse_args()


# create the gas and set initial quiescent state
gas = ct.Solution(init_state.mech)
T1 = init_state.te;
D1 = init_state.den;
gas.TDX = T1, D1, init_state.comp


# Create a reactor called "combustor" 
combustor = ct.IdealGasReactor(gas)
sim = ct.ReactorNet([combustor])



# store intermediate state of the combustor in a CSV file
csvfile='combustor.csv'
outfile = open(csvfile,'w')
csvwriter=csv.writer(outfile, lineterminator="\n")
# the header of the output file
csvwriter.writerow(['time (s)','T (K)','dT/dt (K/s)', 'P (Pa)','rho'] + gas.species_names )

# evolve the combustor by one step and record its state
# this is done with the sim.step(tfinal) which outputs the state after a solver timeste,
# which can be smaller than the overall target finaltime.
# note that sim.advance(tfinal) advaces the reactor to finaltime
gas.TDX = T1, D1, init_state.comp
combustor = ct.IdealGasReactor(gas)
sim = ct.ReactorNet([combustor])
sim.rtol = 1E-5
sim.atol = 1E-12
# settings for SanDiego mech and acetylene combustion
#sim.set_max_time_step(1E-12)
#tfinal = 1E-6
#sim.set_max_time_step(1E-8)
#tfinal = 0.01
tfinal = 1.0
tnow = 0.0
Told = T1
told = 0.0
Tdotold = 0.0
tignition = 0.0
texplosion = 0.0
temp = []
while tnow < tfinal:
    tnow = sim.step()
    temp.append = combustor.thermo.T
# 	if (tnow > told):
# 	    	Tdot = (combustor.T - Told) / (tnow - told)
# 	else:
# 	        Tdot = 0.0
# 	Told = combustor.T
# 	told = tnow
# 	if Tdot > Tdotold:
# 	       Tdotold = Tdot
# 	       tignition = tnow
# 	       texplosion = combustor.T / Tdot
# 	csvwriter.writerow([tnow,  combustor.T, Tdot, combustor.thermo.P,combustor.thermo.density] + list(combustor.thermo.X))
# #print 'output written to '+csvfile
# print 'Ignition time for T1=', T1, ' is ', tignition, ' seconds'
# #print 't_explosion is', texplosion, 'seconds'
# outfile.close()

import matplotlib.pyplot as plt
plt.figure()
plt.plot(temp)
