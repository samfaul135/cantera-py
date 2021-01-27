from SDToolbox import *
from optparse import OptionParser
parser = OptionParser()
OneAtm = 101300	# Pa
GasConstant = 8314
parser.add_option("--pg", type="float", dest="pg", default=OneAtm)
parser.add_option("--te", type="float", dest="te", default=300.0)
parser.add_option("--x", type="string", dest="comp", default='CH4:1,O2:2')
parser.add_option("--mech", type="string", dest="mech", default='gri30_highT.cti')
(init_state, args) = parser.parse_args()


# set initial quiescent state
P1 = init_state.pg;
T1 = init_state.te;
#q = 'CH4:1 O2:2'
#mech = 'gri30_highT.cti'

# compute initial sound speed
gas = Solution(init_state.mech)
gas.TPX = T1, P1, init_state.comp
rho = gas.density
M = gas.mean_molecular_weight
cp0 = 0.001*gas.cp_mass
cv0 = 0.001*gas.cv_mass
gamma = cp0/cv0
c = math.sqrt(gamma*GasConstant*gas.T/M)
print ' '
print 'UNDISTURBED STATE'
print 'Pressure = %.2f Pa' % (P1)
print 'Temperature = %.2f K' % (T1)
print 'Density = %.2f kg/m^3' % (rho)
print 'Sound speed = %.2f m/s' % c
print 'gamma = %.2f' % gamma
print ' '

# compute CJ Speed and CJ/VN states
[cj_speed,R2] = CJspeed(P1, T1, init_state.comp, init_state.mech, 0);   
gas_cj = PostShock_eq(cj_speed, P1, T1, init_state.comp, init_state.mech)
gas_vn = PostShock_fr(cj_speed, P1, T1, init_state.comp, init_state.mech)
Dcj = cj_speed/c
g_cj = gas_cj.cp_mass/gas_cj.cv_mass
g_vn = gas_vn.cp_mass/gas_vn.cv_mass
u_cj = cj_speed - (rho*cj_speed/gas_cj.density)
u_vn = cj_speed - (rho*cj_speed/gas_vn.density)

# get 70% state to find T and P
U_ps = 0.7*cj_speed	# shock speed
gas_ps = PostShock_fr(U_ps, P1, T1, init_state.comp, init_state.mech)
Ps = gas_ps.P/OneAtm
Ts = gas_ps.T
g_ps = gas_ps.cp_mass/gas_ps.cv_mass


print 'CJ STATE'
print 'CJ Speed = %.2f m/s ( Dcj = %.2f )' % (cj_speed,Dcj)
print 'CJ Pressure = %.2f Pa' % (gas_cj.P)
print 'CJ Temperature = %.2f K' % (gas_cj.T)
print 'CJ Density = %.2f kg/m^3' % (gas_cj.density)
print 'CJ Velocity = %.2f m/s' % u_cj
print 'CJ gamma = %.2f' % g_cj
print ' '
print 'VN STATE'
print 'VN Pressure = %.2f Pa' % (gas_vn.P)
print 'VN Temperature = %.2f K' % (gas_vn.T)
print 'VN Density = %.2f kg/m^3' % (gas_vn.density)
print 'VN Velocity = %.2f m/s' % u_vn
print 'VN gamma = %.2f' % g_vn
print ' '
print '70% CJ STATE'
print 'PS Pressure = %.2f atm' % (Ps)
print 'PS Temperature = %.2f K' % (gas_ps.T)
print 'PS Density = %.2f kg/m^3' % (gas_ps.density)
print 'PS gamma = %.2f' % g_ps
print ' '

# Get heat release using VN state gamma
Q = c*c*(Dcj-(1/Dcj))*(Dcj-(1/Dcj))/(2*((g_vn*g_vn) - 1))
print 'Heat release = %.2f J/kg' % Q


# Tune EA and A to give lamda and SL



