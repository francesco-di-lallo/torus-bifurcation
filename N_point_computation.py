# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 23:54:54 2020

@author: Francesco
"""

import matplotlib.pyplot as plt
from ode import Ode
from bmk import Bmk
import numpy as np

def jacobian(PE, ox, oy, step = 1e-5):
    ##PE_x
    bmk = Bmk({'ox':ox+step, 'oy':oy})
    ode = Ode(bmk.ode)
    
    PE_x = ode.pontryagin_energy()
    
    bmk = Bmk({'ox':ox, 'oy':oy+step})
    ode = Ode(bmk.ode)
    
    PE_y = ode.pontryagin_energy()
    
    return np.array([[(PE_x[0]-PE[0])/step, (PE_y[0]-PE[0])/step],
            [(PE_x[1]-PE[1])/step, (PE_y[1]-PE[1])/step]])

par = [-0.00964988,  0.92704329]
progress = [tuple(par)]
PE_progress = []

for i in range(100):
    print(i)
    bmk = Bmk({'ox':par[0], 'oy':par[1]})
    ode = Ode(bmk.ode)
    
    PE = ode.pontryagin_energy()
    
    jac = jacobian(PE, par[0], par[1])
    
    par -= 0.1*np.matmul(np.linalg.inv(jac), np.array(PE))
    progress.append(tuple(par))
    PE_progress.append(sum(PE))
    
#pickle.dump(progress, open('progress_08.p','wb'))   
#plt.scatter(*zip(*progress))
plt.plot(PE_progress)
plt.show()