# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 15:35:23 2020

@author: Francesco
"""

import numpy as np
from bmk import Bmk
from ode import Ode
import matplotlib.pyplot as plt
from itertools import product
import pickle
"""
def f(x,y,ox,oy):
    x_rhs = -ox+np.cos(2*np.pi*y)+0.1*np.cos(2*np.pi*x)
    y_rhs = -oy+np.sin(2*np.pi*y)+0.1*np.sin(2*np.pi*x)
    
    return x_rhs, y_rhs

def euler_heun_iter(x,y,ox,oy,h):
    (x_dot, y_dot) = f(x,y,ox,oy)
        
    (x_tilde, y_tilde) = (x+float(h)*x_dot, y+float(h)*y_dot)
    
    (x_tilde_dot, y_tilde_dot) = f(x_tilde,y_tilde,ox,oy)
    
    (x_next, y_next) = (x+float(h/2)*(x_dot+x_tilde_dot),
                        y+float(h/2)*(y_dot+y_tilde_dot))
    
    return x_next, y_next

def iterate(ox, oy, seed):

    x = 0

    curve_x = [x]
    curve_y = [seed]
    
    (x,y) = euler_heun_iter(x,seed,ox,oy,1e-2)

    curve_x.append(x)
    curve_y.append(y)
    looper = 0
    while x>-1 and looper < 500:# curve_x[-1]<curve_x[-2]:
        looper += 1
        (x,y) = euler_heun_iter(x,y,ox,oy,1e-2)
        curve_x.append(x)
        curve_y.append(y)
    if abs(x+1) < 1e-1:    
        return y
    else:
        #delta_H.append(None)
        return None  


def exist_po(ox, oy):
    y = 0.5
    for _ in range(20):
        y = iterate(ox, oy, y%1)
        #print(y)
        if y == None:
            break
        #print(y)
    
    if y != None:
        return True
        #(np.sign(iterate(ox, oy, y-1e-3)-y))
    else:
        return False

y = 0
#ox = -0.125
ox = 0.2
oy = 0.9

for oy in np.linspace(0.95,0.98,20):
    
    print(exist_po(ox, oy), (ox, oy))

rot_snp_ox = [ox]
rot_snp_oy = [oy]
mesh = 1e-2

while ox < 0.55:
    
    oy = rot_snp_oy[-1]
    
    print(ox, oy)
    
    po_bool = exist_po(ox, oy)
    if po_bool == True:
        step = mesh
    else:
        step = -mesh
    while po_bool == exist_po(ox, oy):
        print("Testing ",oy)
        oy += step
    ox += 1e-2
    rot_snp_ox.append(ox)
    rot_snp_oy.append(oy)
"""

rot_snp_ox = list(np.linspace(-0.1,0.55,50))
rot_snp_oy = [
    0.9054,
    0.9087,
    0.9120,
    0.9152,
    0.9184,
    0.9215,
    0.9246,
    0.9276,
    0.9306,
    0.9336,
    0.9363,
    0.9391,
    0.9419,
    0.9445,
    0.9471,
    0.9497,
    0.9522,
    0.9546,
    0.9569,
    0.9592,
    0.9614,
    0.9635,
    0.9655,
    0.9675,
    0.9694,
    0.9712,
    0.9729,
    0.9746,
    0.9762,
    0.9777,
    0.9791,
    0.9806,
    0.9818,
    0.9830,
    0.9842,
    0.9854,
    0.9865,
    0.9874,
    0.9881,
    0.9890,
    0.9898,
    0.9907,
    0.9912,
    0.9919,
    0.9927,
    0.9934,
    0.9941,
    0.9947,
    0.9953,
    0.9959
    ]
pickle.dump(rot_snp_ox, open('rot_snp_plus_ox.p', 'wb'))
pickle.dump(rot_snp_oy, open('rot_snp_plus_oy.p', 'wb'))
