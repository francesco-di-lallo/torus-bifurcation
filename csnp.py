# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 13:46:02 2020

@author: Francesco
"""


from bmk import Bmk
from ode import Ode
import numpy as np
import matplotlib.pyplot as plt
import pickle

def transform(x, y, c, s, centre):
    x -= centre[0]['x']
    y -= centre[0]['y']
    return (c*x+s*y, c*y-s*x)

def periodic_counter(ox, oy, granularity = 200):
    bmk = Bmk({'ox':ox, 'oy':oy})
    ode = Ode(bmk.ode)
    
    centre = ode.get_centre()
    saddle = ode.get_saddle_point()
    
    if centre[0]['x'] > saddle[0]['x']:
    
        theta = np.arctan((centre[0]['y']-saddle[0]['y'])/
                      (centre[0]['x']-saddle[0]['x']))
        
    else:
        
        theta = -np.pi + np.arctan((centre[0]['y']-saddle[0]['y'])/
                      (centre[0]['x']-saddle[0]['x']))
    
    c = np.cos(theta)
    s = np.sin(theta)
    
    granularity = 80
    
    
    H_axis_x = np.linspace(centre[0]['x'], saddle[0]['x'], granularity)
    H_axis_y = np.linspace(centre[0]['y'], saddle[0]['y'], granularity)
    
    H_axis = list(zip(H_axis_x, H_axis_y))
    
    delta_H = []
    
    
    for (x_0, y_0) in H_axis:
        
        pt_dict = {'x':x_0, 'y':y_0}
        (x_0, y_0) = transform(x_0, y_0, c, s, centre)
        curve_x = [x_0]
        curve_y = [y_0]
        not_crossed = True
        looper = 0
        while not_crossed and looper < 50000:
            looper += 1
            pt_dict = ode._euler_iter(pt_dict, h=5e-3)
            pt_x, pt_y = transform(pt_dict['x'], pt_dict['y'], c, s, centre)
            
            if curve_y[-1] < -1e-12 and pt_y > 1e-12:
                not_crossed = False
            
            curve_x.append(pt_x)
            curve_y.append(pt_y)
        if np.abs(curve_y[-1])<1e-1 and curve_x[-1]<1e-2:
            delta_H.append(x_0-pt_x)
        else:
            #print('Escaped bounded region')
            break
       # plt.plot(curve_x, curve_y)
    #plt.show()
    po_count = 0
    index = []
    for H in range(10,len(delta_H)-1):
        if np.sign(delta_H[H]) != np.sign(delta_H[H-1]):
            po_count += 1
            index.append(H)

         
    if po_count == 2:
        
        return po_count, delta_H
    
    else:
        
        return po_count, delta_H
    
"""
ox = np.linspace(0.035, 0.04,10)
oy = [0.95]*10

po = []

looper = 0
for (x, y) in zip(ox, oy):
    looper += 1
    print(looper, '(ox, oy)', (round(x,6), round(y,6)))
    p, d = periodic_counter(x,y)
    
    po.append(p)
    
    print("Periodic Orbits Detected: ",p)
    
    plt.plot(d)
        

plt.show()
"""
csnp = pickle.load(open('two_po.p', 'rb'))
two_x = [po[0] for po in csnp]
two_y = [po[1] for po in csnp]

csnp_x = []
csnp_y = []

for oy in set(two_y):
    poss = []
    for pair in csnp:
        if pair[1] == oy:
            poss.append(pair[0])
    csnp_y.append(oy)
    csnp_x.append(max(poss))
"""
(ox, oy) = (0.039, 0.95)
two_po = [(ox, oy)]
visited = [(ox, oy)]

mesh = 1e-4
sig_fig = 4

queue = [(ox, oy-mesh), (ox-mesh, oy), (ox+mesh, oy), (ox, oy+mesh)]

while len(queue) != 0:
    (ox, oy) = queue.pop(-1)
    ox = round(ox,sig_fig)
    oy = round(oy,sig_fig)
    if (ox, oy) not in visited:
        p, dist = periodic_counter(ox,oy)
        
        if p==2:
            two_po.append((ox, oy))
            queue += [(ox, round(oy-mesh,sig_fig)), (round(ox-mesh,sig_fig), oy),
                      (round(ox+mesh,sig_fig), oy), (ox, round(oy+mesh,sig_fig))]
    print(len(queue), ox, oy)
    visited.append((ox, oy))
"""
"""
D points are at :
    (ox, oy) = (0.00932643537002703, 0.9319331314841488)
    
    (ox, oy) = (0.09999371598888401, 0.9999937161665724)
    
csnp : (ox,oy) \contains (0.04169590947368421, 0.95)
"""