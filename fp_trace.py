# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 21:57:51 2020

@author: Francesco
"""


import matplotlib.pyplot as plt
from ode import Ode
from bmk import Bmk
import numpy as np

ox = np.linspace(-1e-02, 0, 100)
tr_saddle = []
tr_node = []
#-0.009649961793522382], [0.9270432434144844
for i in ox:
    bmk = Bmk({'ox': i, 'oy': 0.9270432434144844})
    ode = Ode(bmk.ode)
    fp = ode.get_fixed_points()
    fp_type = ode.fp_eigen(fp[0])[0]
    if fp_type == 'saddle point':
        print('case saddle {}'.format(i))
        tr_saddle.append(ode.jacobian(fp[0])[0][0]+ode.jacobian(fp[0])[1][1])
        tr_node.append(ode.jacobian(fp[1])[0][0]+ode.jacobian(fp[1])[1][1])
    else:
        print('case node {}'.format(i))
        tr_saddle.append(ode.jacobian(fp[1])[0][0]+ode.jacobian(fp[1])[1][1])
        tr_node.append(ode.jacobian(fp[0])[0][0]+ode.jacobian(fp[0])[1][1])
        
#figure, axes = plt.subplots(nrows=2, ncols=1)        
plt.plot(ox, tr_saddle)
#plt.plot(ox, [0]*len(ox))
plt.plot([-0.009649961793522382]*2,[min(min(tr_saddle),min(tr_node)),
                          max(max(tr_saddle),max(tr_node))])
plt.plot(ox, [0]*len(ox), c='k')
plt.plot(ox, tr_node)
plt.show()