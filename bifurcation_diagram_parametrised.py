# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 17:55:37 2020

@author: Francesco
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
from bmk import Bmk
from ode import Ode

eps = 1e-1

bound = np.arctan(eps-np.sqrt(eps**2+1))/np.pi
t_tr0_centre = np.linspace(bound,-bound,1000)
t_tr0_saddle = np.linspace(-bound, 1-bound,1000)
t_sne = np.linspace(0,2*np.pi, 1000)

def inner_sne(t, eps):
    ox = (1-eps)*np.cos(t)
    oy = (1-eps)*np.sin(t)
    
    return ox, oy

def outer_sne(t, eps):
    ox = (1+eps)*np.cos(t)
    oy = (1+eps)*np.sin(t)
    
    return ox, oy

def trace_zero_upper(t, eps):
    ox = eps*np.sin(2*np.pi*t) + eps*np.cos(2*np.pi*t)
    oy = np.sqrt(1-(eps*np.sin(2*np.pi*t))**2) + eps*np.sin(2*np.pi*t)
    
    return ox, oy

def trace_zero_lower(t, eps):
    ox = eps*np.sin(2*np.pi*t) + eps*np.cos(2*np.pi*t)
    oy = -1*np.sqrt(1-(eps*np.sin(2*np.pi*t))**2) + eps*np.sin(2*np.pi*t)
    
    return ox, oy


def scale(ox, oy, eps):
    if len(ox) != len(oy):
        print("Omega_x and Omega_y do not have the same length")
        return None

    rho = [(1/eps)*(norm(x,y)-1) for (x,y) in zip(ox, oy)]
    alphat = [(1/eps)*(np.pi/2 - np.arccos(x)/norm(x,y)) for (x,y) in zip(ox, oy)]
    
    return alphat, rho

def norm(a,b):
    return np.sqrt(a**2 + b**2)


i_x, i_y = inner_sne(t_sne, eps)
o_x, o_y = outer_sne(t_sne, eps)
U_centre_x, U_centre_y = trace_zero_upper(t_tr0_centre, eps)
U_saddle_x, U_saddle_y = trace_zero_upper(t_tr0_saddle, eps)
#L_x, L_y = trace_zero_lower(t_tr0, eps)


"""
plt.plot(np.linspace(0,1,len(t_tr0)), ly)
plt.title("First Lyapunov Coefficient of Curve of Centre")
plt.xlabel("Arclength parameter of trace zero")
plt.show()
"""
"""
rhc_oy = pickle.load(open('rhc_oy.p', 'rb'))
rhc_ox = pickle.load(open('rhc_ox.p', 'rb'))
"""
rhc_oy = pickle.load(open('rhc_oy_heun.p', 'rb'))
rhc_ox = pickle.load(open('rhc_ox_heun.p', 'rb'))

rhc_oy_lift = pickle.load(open('rhc_oy_lift_heun.p', 'rb'))
rhc_ox_lift = pickle.load(open('rhc_ox_lift_heun.p', 'rb'))

chc_ox = pickle.load(open('chc_ox_heun.p', 'rb'))
chc_oy = pickle.load(open('chc_oy_heun.p', 'rb'))


csnp_x = pickle.load(open('csnp_x.p', 'rb'))
csnp_y = pickle.load(open('csnp_y.p', 'rb'))
csnp_x
plt.plot(csnp_x, csnp_y, color='#ff8000', zorder=4, linewidth=1.2)

rot_snp_ox = pickle.load(open('rot_snp_ox.p', 'rb'))
rot_snp_oy = pickle.load(open('rot_snp_oy.p', 'rb'))
rot_snp_plus_ox = pickle.load(open('rot_snp_plus_ox.p', 'rb'))
rot_snp_plus_oy = pickle.load(open('rot_snp_plus_oy.p', 'rb'))

rot_snp_ox += [-0.6]
rot_snp_oy += [0.995]



rot_snp_plus_ox += [0.6]
rot_snp_plus_oy += [0.995]

rot_snp_plus_oy = [oy + 0.0006 for oy in rot_snp_plus_oy]
rot_snp_plus_ox.insert(0,-0.124)
rot_snp_plus_oy.insert(0,0.8994)
rot_snp_oy[0] = 0.9517
rot_snp_ox[0] = -0.1355
plt.plot(rot_snp_plus_ox, rot_snp_plus_oy, c='g')
plt.plot(rot_snp_ox, rot_snp_oy, c='g')

"""
i_alphat, i_rho = scale(i_x, i_y, eps)
o_alphat, o_rho = scale(o_x, o_y, eps)
U_centre_alphat, U_centre_rho = scale(U_centre_x, U_centre_y, eps)
U_saddle_alphat, U_saddle_rho = scale(U_saddle_x, U_saddle_y, eps)
#L_alphat, L_rho = scale(L_x, L_y, eps)
rhc_alphat, rhc_rho = scale(rhc_ox, rhc_oy, eps)
rhc_lift_alphat, rhc_lift_rho = scale(rhc_ox_lift, rhc_oy_lift, eps)
"""


"""
necklace_alphat, necklace_rho = scale([-0.009649961793522382],
                                      [0.9270432434144844], 0.1)

plt.plot(i_alphat, i_rho, c='b')
plt.plot(o_alphat, o_rho, c='b')
plt.plot(U_centre_alphat, U_centre_rho, c='k')
plt.plot(U_saddle_alphat, U_saddle_rho, c='k')
#plt.plot(L_alphat, L_rho, c='k')
plt.plot(rhc_alphat, rhc_rho, c='r')
plt.plot(rhc_lift_alphat, rhc_lift_rho, c='r')


plt.scatter([U_centre_alphat[0], U_centre_alphat[-1]],
            [U_centre_rho[0], U_centre_rho[-1]], marker = 's',zorder=3)

plt.scatter(necklace_alphat, necklace_rho, c='k', marker = 's',zorder=3)

plt.xlabel('alpha tilde')
plt.ylabel('rho')
plt.xlim(-3.2,2.8)
plt.ylim(-1.1,1.1)

plt.savefig('Bifurcation Diagram Scale Par.png', dpi = 2000)
plt.show()
"""
PT_SIZE = 6

plt.plot(i_x, i_y, c='b', linewidth=1)
plt.plot(o_x, o_y, c='b', linewidth=1)
plt.plot(U_centre_x, U_centre_y, c='k', linewidth=1)
plt.plot(U_saddle_x, U_saddle_y, c='k', linestyle='dashed', linewidth=1)
#plt.plot(L_x, L_y, c='k')
plt.plot(rhc_ox, rhc_oy, c='r', linewidth=1)
buffer = 150
plt.plot(rhc_ox_lift[:buffer], rhc_oy_lift[:buffer], c='r', linewidth=1)
chc_ox = [ox-4e-4 for ox in chc_ox]
plt.plot(chc_ox, chc_oy, c='r', linestyle='dashed', linewidth=1)


plt.scatter([U_centre_x[0], U_centre_x[-1]],
            [U_centre_y[0], U_centre_y[-1]],
            s=PT_SIZE, c='k',zorder=3)


plt.scatter([-0.0103], [0.9270432434144844],
            s=PT_SIZE,c='k',zorder=3)

plt.scatter([rhc_ox[0], rhc_ox[-1], rhc_ox_lift[0], rhc_ox_lift[buffer-1]],
            [rhc_oy[0], rhc_oy[-1], rhc_oy_lift[0], rhc_oy_lift[buffer-1]],
            s=PT_SIZE,c='k', zorder=3)

plt.scatter(0.0174, 0.93664,
            s=PT_SIZE, c='k', zorder=5)

plt.scatter(csnp_x[-1], csnp_y[-1],
            s=PT_SIZE, c='k', zorder=5)

plt.scatter(-0.124, 0.8994,
            s=PT_SIZE, c='k', zorder=5)
plt.scatter(-0.0113, 0.92727,
            s=PT_SIZE, c='k', zorder=5)
plt.scatter(rot_snp_ox[0], rot_snp_oy[0],
            s=PT_SIZE, c='k', zorder=5)

plt.annotate('$N$', (-0.009649961793522382, 0.9270432434144844),
             xytext=(-0.017, 0.938))

plt.annotate('$H$', (-0.009649961793522382, 0.9270432434144844),
             xytext=(-0.042, 0.938))

plt.annotate('$B^+$', (U_centre_x[0], U_centre_y[0]),
             xytext=(-0.1,0.85)) ###
plt.annotate('$B^-$', (U_centre_x[-1], U_centre_y[-1]),
             xytext=(0.1,1.11)) ###

plt.annotate('$Z_a^+$',(rhc_ox[0], rhc_oy[0]),
             xytext=(-0.28,0.884))

plt.annotate('$Z_a^+$',(rhc_ox[-1], rhc_oy[-1]),
             xytext=(0.53,0.898)) ###

plt.annotate('$Z_b^-$',(rhc_ox_lift[0], rhc_oy_lift[0]),
             xytext=(-0.53, 0.925)) ###

plt.annotate('$Z_b^-$',(rhc_ox_lift[buffer], rhc_oy_lift[buffer]),
             xytext=(0.3,0.865))###

plt.annotate('$K$',(rhc_ox_lift[buffer], rhc_oy_lift[buffer]),
             xytext=(-0.17,0.97))###

plt.annotate('$K$',(rhc_ox_lift[buffer], rhc_oy_lift[buffer]),
             xytext=(-0.17,0.9))###

plt.annotate('$D$', (0.015,0.94),
             xytext=(0.015,0.9489))###

plt.annotate('$D$', (rot_snp_plus_ox[0], rot_snp_plus_oy[0]),
             xytext=(0.09,0.97))###

#plt.scatter(0.041703706666666666, 0.95)
plt.xlabel('$\Omega_x$')
plt.ylabel('$\Omega_y$')

plt.axis('equal')

plt.xlim(-0.6,0.6)

plt.ylim(0.7,1.11)
"""
plt.xlim(-0.02, 0.03)
plt.ylim(0.92,0.945)
"""
plt.savefig('Bifurcation Diagram all.png', dpi = 2000)
plt.show()
