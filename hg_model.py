from __future__ import print_function
import numpy as np
from boxpy import *
import pylab as pl
from inifile import ini
from fromfile import froms


compartments = ['atm','land_arm','land_fast','land_slow',
                'nor_at','deep_nor_at','sfc_at','int_at',
                'deep_at','bot_at','med_sea','deep_med','sou_oc',
                'deep_sou_oc','int_pac','deep_pac','sfc_pac',
                'nor_pac']

K = construct_K(froms,compartments)

geogen = np.zeros(len(compartments))
geogen[0] = 90.

def s(tt):
    """Emissions forcing function"""
    extra = np.zeros_like(geogen)
    if tt < 1050:
        pass
    elif tt < 1850:
        extra[0] += 500.
    elif tt < 1950:     
        extra[0] += 500 + 2000*np.exp((-(tt-1890)**2)/(20**2))
    elif tt < 2000:
        extra[0] += 500. + 1000*(tt-1950)/50.
    elif tt < 2100:
        extra[0] += 1500. + 3000*(tt-2000)/50.

    return geogen + extra

c0 = np.array([ini[name] for name in compartments])

t = np.linspace(-2000,2050,4051)
c = solve_odeint(get_f_ddt(K,s),c0,t)


cc=1
pl.subplot(2,2,cc)
pl.plot(t,[s(tt)[0] for tt in t], label='emis')
pl.legend()
pl.xlim(1840,2050)

cc+=1
for i in [0,3,-1]:
    pl.subplot(2,2,cc)
    pl.plot(t[:],c[:,i],label=compartments[i])
    pl.xlim(1800,2050)
    pl.legend()
    cc += 1


cs = np.sum(c,axis=1)
#pl.figure()
#pl.plot(t,cs)
print((cs[-1]-cs[0])/4050.)

pl.show()
