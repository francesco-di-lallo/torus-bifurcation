# -*- coding: utf-8 -*-
"""
Created on Fri May 29 18:43:12 2020

@author: Francesco
"""


import matplotlib.pyplot as plt
from ode import Ode
from bmk import Bmk
import numpy as np

bmk = Bmk({'ox': -0.009, 'oy': 0.927})
ode = Ode(bmk.ode)

fp = ode.get_fixed_points()
if ode.fp_eigen(fp[0]) in ['hyperbolic repeller', 'hyperbolic attractor']:
    centre = fp[0]
    saddle = fp[1]
else:
    centre = fp[1]
    saddle = fp[0]

H_axis_x = np.linspace(centre['x'],saddle['x']+1, 20)
H_axis_y = np.linspace(centre['y'],centre['y'], 20)

H_axis = list(zip(H_axis_x, H_axis_y))

