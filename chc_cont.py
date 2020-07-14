# -*- coding: utf-8 -*-
"""
Created on Thu May 21 23:10:25 2020

@author: Francesco
"""

from ode import Ode
from bmk import Bmk
import numpy as np
import matplotlib.pyplot as plt
import pickle

def learning(ox_seed, oy, direction):
    ox = ox_seed 
    learning_rate = 5e-1
    estimates = [ox_seed]
    ox -= learning_rate*newton_raphson(ox, oy, direction)
    estimates.append(ox)
    i = 0
    while np.abs(estimates[-1]-estimates[-2])>5e-6 and i<50:
        print('Estimate:',round(ox,6),
              ' delta:', round(estimates[-1]-estimates[-2],6))
        ox -= learning_rate*newton_raphson(ox, oy, direction)
        estimates.append(ox)
        i+=1
    return estimates[-1]

def newton_raphson(ox, oy, direction, step = 1e-8):
    
    bmk = Bmk({'ox':ox, 'oy':oy})
    ode = Ode(bmk.ode)
    
    PE = ode.homoclinic_dist(direction, 'contractible')
    
    
    bmk_incr = Bmk({'ox':ox+step, 'oy':oy})
    ode_incr = Ode(bmk_incr.ode)
    
    PE_incr = ode_incr.homoclinic_dist(direction, 'contractible')
    return PE/((PE_incr-PE)/step)


N = [-0.010115, 0.9270432434144844]
def propagate(N):
    chc_ox = [N[0]]
    chc_oy = [N[1]]
    oy = N[1]
    for i in np.arange(oy,1.1,0.0005):
        print(i)
        ox_seed = chc_ox[-1]
        oy_seed = i
        try:
            chc_ox.append(learning(ox_seed, oy_seed, -1))
            chc_oy.append(i)
        except:
            print("done goofed")
    
    for i in np.arange(oy,0.894,-0.0005):
        print(i)
        ox_seed = chc_ox[0]
        oy_seed = i
        try:
            chc_ox.insert(0, learning(ox_seed, oy_seed,1))
            chc_oy.insert(0, i)
        except:
            print("done goofed")
    
    plt.plot(chc_ox, chc_oy)
    pickle.dump(chc_ox, open('chc_ox_heun.p', 'wb'))
    pickle.dump(chc_oy, open('chc_oy_heun.p', 'wb'))
    return chc_ox, chc_oy
