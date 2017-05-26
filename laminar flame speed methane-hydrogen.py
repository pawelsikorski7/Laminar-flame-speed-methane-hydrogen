"""
A freely-propagating, premixed flat flame with multicomponent
transport properties.
"""
import cantera as ct
import numpy as np
import csv
import os

#print ct.__file__

p = ct.one_atm
tburner = 300
gas = ct.Solution('gri30.cti')
gas.transport_model = 'Mix'
mdot = 0.6

initial_grid = np.linspace(-0.01, 0.4, 9)

tol_ss = [1.0e-5, 1.0e-10]  # [rtol atol] for steady-state problem 
tol_ts = [1.0e-5, 1.0e-10]  # [rtol atol] for time stepping loglevel = 1 # amount of diagnostic output (0 to 8)
refine_grid = True  # 'True' to enable refinement, 'False' to disable

# Mole Fractions
EQUI = 1.4
EQUI_END = 2.7
phi_step = 0.1
phi_steps = (np.ceil((EQUI_END - EQUI)/phi_step))+1
n=0
phi=[0]*phi_steps
fs=[0]*phi_steps

io2 = gas.species_index('O2'); #Index of O2 in mix
in2 = gas.species_index('N2'); #Index of N2 in mix
ih2 = gas.species_index('H2'); #Index of N2 in mix
ich4 = gas.species_index('CH4'); #Index of N2 in mix
iar  = gas.species_index('AR'); #Index of AR in mix 

for EQUI in np.arange(EQUI, EQUI_END+0.001, phi_step):
        F1_Frac = 1.0;#Fraction of the first fuel in the blend
        F2_Frac = 0.0;#Fraction of the second fuel in the blend
        F3_Frac = 0.0;#Fraction of the third fuel in the blend
        #c = 2; #Number of moles of O2 for CH4
        #c= 1.25 #Number of moles for 50-50
        c=0.5 #Number of moles for H2

        comp = [0]*gas.n_species;
        comp[ich4] = F2_Frac * EQUI
        comp[ih2] = F1_Frac * EQUI
        comp[io2] = c
        comp[in2] = 3.76*c * 0.78/0.79  
        comp[iar] = 3.76*c * 0.01/0.79  

        gas.TPX = tburner, p, comp

        # Flame object
        f = ct.FreeFlame(gas, initial_grid)
        f.flame.set_steady_tolerances(default=tol_ss)
        f.flame.set_transient_tolerances(default=tol_ts)
        f.set_refine_criteria(ratio=2.7, slope=0.06, curve=0.12)
        f.set_grid_min(1e-9)

        # Set properties of the upstream fuel-air mixture
        f.inlet.T = tburner
        f.inlet.X = comp
        f.inlet.mdot = mdot

        f.transport_model = 'Mix'
        f.set_time_step(1e-5, [2, 5, 10])
        f.energy_enabled = False
        try:
                f.solve(loglevel=loglevel, refine_grid=False)
        except Exception:
                print ("failed initial solve at phi = ", EQUI)
        f.energy_enabled = True

        try:
                f.solve(loglevel=loglevel, refine_grid=refine_grid)
        except Exception:
                print ("failed mix solve at phi = ", EQUI)
        print ("\n***************MIX ENERGY {0:7f} {1:7f}***************".format(EQUI, f.u[0]))

        f.transport_model = 'Multi'
        try:
                f.solve(loglevel=loglevel, refine_grid=refine_grid)
        except Exception:
                print ("failed multi solve at phi = ", EQUI)
        print ("\n***************MULTI ENERGY {0:7f} {1:7f}***************".format(EQUI, f.u[0]))

        f.soret_enabled = True
        try:
                f.solve(loglevel=loglevel, refine_grid=refine_grid)
        except Exception:
                print ("failed soret solve at phi = ", EQUI)
        f.write_csv('profile-phi{:.1f}.csv'.format(EQUI), quiet=False)
        print ("\n***************MULTI ENERGY SORET {0:7f} {1:7f}***************".format(EQUI, f.u[0]))

        phi[n]=EQUI
        fs[n] = f.u[0]
        n=n+1
