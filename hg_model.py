from __future__ import print_function
import numpy as np
from boxpy import *
import pylab as pl

compartments = ['atm','land_arm','land_fast','land_slow',
                'nor_at','deep_nor_at','sfc_at','int_at',
                'deep_at','bot_at','med_sea','deep_med','sou_oc',
                'deep_sou_oc','int_pac','deep_pac','sfc_pac',
                'nor_pac']


froms = {
    'atm':{
        'atm': 0.,
        'land_arm': 0.05,
        'land_fast': 0.4,
        'land_slow': 0.08,
        'nor_at': 0.7*0.037,
        'deep_nor_at': 0.,
        'sfc_at': 0.7*0.17,
        'int_at': 0.7*0.072,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.7*0.0053,
        'deep_med': 0.,
        'sou_oc': 0.7*0.027,
        'deep_sou_oc': 0.,
        'int_pac': 0.7*0.176,
        'deep_pac': 0.,
        'sfc_pac': 0.7*0.42,
        'nor_pac': 0.7*0.09,
        },
    
    'land_arm':{
        'atm': 2.e-4,
        'land_arm': 0.,
        'land_fast': 8.e-5,
        'land_slow': 0.,
        'nor_at': 3.e-5*0.02,
        'deep_nor_at': 0.,
        'sfc_at': 3.e-5*0.36,
        'int_at': 3.e-5*0.02,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 3.e-5*0.04,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 3.e-5*0.26,
        'nor_pac': 3.e-5*0.30,
        },
    
    'land_fast':{
        'atm': 0.2,
        'land_arm': 1.e-3,
        'land_fast': 0.,
        'land_slow': 0.03,
        'nor_at': 0.04*0.02,
        'deep_nor_at': 0.,
        'sfc_at': 0.04*0.36,
        'int_at': 0.04*0.02,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.04*0.04,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.04*0.26,
        'nor_pac': 0.04*0.30,
        },
    
    'land_slow':{
        'atm': 7.e-3,
        'land_arm': 1.e-5,
        'land_fast': 6.e-3,
        'land_slow': 0.,
        'nor_at': 3.e-4*0.02,
        'deep_nor_at': 0.,
        'sfc_at': 3.e-4*0.36,
        'int_at': 3.e-4*0.02,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 3.e-4*0.04,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 3.e-4*0.26,
        'nor_pac': 3.e-4*0.30,
        },
    
    'nor_at':{
        'atm': 6.e-3,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.02,
        'sfc_at': 0.,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'deep_nor_at':{
        'atm': 0.,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': -3.e-4,
        'sfc_at': 0.,
        'int_at': 0.,
        'deep_at': 0.02,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'sfc_at':{
        'atm': 0.1,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.03,
        'deep_nor_at': 0.,
        'sfc_at': 0.,
        'int_at': 0.1,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 1.e-3,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'int_at':{
        'atm': 8.e-3,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 2.e-3,
        'int_at': 0.,
        'deep_at': 8.e-3,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'deep_at':{
        'atm': 0.,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 0.,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': 3.e-3,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 5.e-3,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'bot_at':{
        'atm': 0.,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 4.e-3,
        'sfc_at': 0.,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': -7.e-3,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'med_sea':{
        'atm': 0.2,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 0.,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.05,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'deep_med':{
        'atm': 0.,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 0.,
        'int_at': 8.e-3,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': -1.e-3,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'sou_oc':{
        'atm': 1.e-3,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 0.02,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.02,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'deep_sou_oc':{
        'atm': 0.,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 0.,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': 6.e-3,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.01,
        'deep_sou_oc': -2.e-4,
        'int_pac': 0.,
        'deep_pac': 0.02,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'int_pac':{
        'atm': 0.01,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 0.,
        'int_at': 7.e-4,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 4.e-3,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'deep_pac':{
        'atm': 0.,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 0.,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 3.e-4,
        'deep_pac': -9.e-4,
        'sfc_pac': 0.,
        'nor_pac': 7.e-4,
        },
    
    'sfc_pac':{
        'atm': 0.1,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 0.,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 9.e-3,
        'deep_sou_oc': 0.,
        'int_pac': 0.07,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'nor_pac':{
        'atm': 0.02,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 9.e-3,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 4.e-3,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    }


    
ini = {'atm': 487,
       'land_arm': 81000,
       'land_fast': 850,
       'land_slow': 4850,
       'nor_at': 2409,
       'deep_nor_at': 3550,
       'sfc_at': 526,
       'int_at': 6417.,
       'deep_at': 15260.,
       'bot_at': 5900.,
       'med_sea': 190.,
       'deep_med': 1500.,
       'sou_oc': 2000.,
       'deep_sou_oc': 3020.,
       'int_pac': 5400.,
       'deep_pac': 50810.,
       'sfc_pac': 260.,
       'nor_pac': 2500.,
        }

#K = np.zeros((len(compartments),len(compartments)))
#for i,namei in enumerate(compartments):
#    for j,namej in enumerate(compartments):
#        K[j,i] += froms[namei][namej]
#        if i != j:
#            K[i,i] -= froms[namei][namej]

K = construct_K(froms,compartments)

print(np.max(K))

geogen = np.zeros(len(compartments))
geogen[0] = 90.




def s(tt):
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
