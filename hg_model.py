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
                'deep_at','bot_at','med_sea','deep_med','sou_oc',
                'deep_sou_oc','int_pac','deep_pac','sfc_pac',
                'nor_pac']

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
pl.xlim(1840,2050)

cc+=1
for i in [0,3,-1]:
    pl.subplot(2,2,cc)
    pl.plot(t[:],c[:,i],label=compartments[i])
    pl.xlim(1800,2050)
    pl.legend(loc='lower right')
    cc += 1


cs = np.sum(c,axis=1)

pl.show()
