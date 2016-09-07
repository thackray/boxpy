"""Example use of boxpy, mercury biogeochem example."""

from __future__ import print_function
import numpy as np
from boxpy import *
import pylab as pl
from inifile import ini
from fromfile import froms
from scipy.interpolate import interp1d


# Name the compartments
compartments = ['atm','land_arm','land_fast','land_slow',
                'nor_at','deep_nor_at','sfc_at','int_at',
                'deep_at','bot_at','med_sea','deep_med',
                'sou_oc','deep_sou_oc','int_pac','deep_pac',
                'sfc_pac','nor_pac']

comp_volumes = [1.,1.,1.,1.,
                2.9e16, 3.6e16, 1.8e16, 1.e17, 
                1.5e17, 4.9e16, 7.5e14, 3.0e15,
                1.7e16, 3.5e16, 2.5e17, 5.4e17,
                4.4e16, 4.1e16]

comp_surfs = [1., 1., 1., 1.,
              1.96e13, 1.96e13, 6.16e13, 2.04e13,
              8.2e13, 8.2e13, 2.5e12, 2.5e12,
              1.15e13, 1.15e13, 5e13, 2.25e14,
              1.5e14, 2.7e13]   


#Box volumes
#2.9e16, 3.6e16, 1.8e16, 1.e17, 
#1.5e17, 4.9e16, 7.5e14, 3.0e15,
#4.1e16, 4.4e16, 2.5e17, 5.4e17, 
#1.7e16, 3.5e16

# Using the dictionary of flows in froms, make flow matrix
K = construct_K(froms,compartments)

# will add geogenic emissions to anthro 
geogen = np.zeros(len(compartments))
geogen[0] = 90.

# Load anthro emissions from a file and combine with geogenic 
#emisfile = np.genfromtxt('AnthroEmissAllTime_20120112_2.txt')
emisfile = np.genfromtxt('anthro_alltime_zero.txt')
eyears = emisfile[:,0]
estreets = np.zeros((len(eyears),len(compartments)))
estreets[:,0] = emisfile[:,1]
estreets += geogen

# Use the function constructor to make the forcing function
s = get_s_from_timeseries(eyears,estreets)

# initial conditions from inifile
c0 = np.array([ini[name] for name in compartments])

# Use the function constructor to make the ddt function
f_ddt = get_f_ddt(K,s)

# Define time axis for output
t = np.linspace(-2000,2050,4051)

# Solve!
c = solve_odeint(f_ddt,c0,t)

# Plot results
pl.figure(figsize=(12,9))
cc=1
pl.subplot(2,2,cc)
pl.plot(t,[s(tt)[0] for tt in t], label='emis')
pl.legend(loc='lower right')
pl.xlim(1800,2050)

cc+=1
for i in [0,3,-1]:
    pl.subplot(2,2,cc)
    pl.plot(t[:],c[:,i],label=compartments[i])
    pl.xlim(1800,2050)
    pl.legend(loc='lower right')
    cc += 1


cs = np.sum(c,axis=1)

# Lifetimes in compartments for 2015:
c2015 = c[2015+2000,:]
lossK = np.diagonal(K)
print(lossK)
dcdt_loss = -1*c2015*lossK
for i,comp in enumerate(compartments):
    print(comp, c2015[i]/dcdt_loss[i])

# Ocean concentrations for 2015:


pl.show()
