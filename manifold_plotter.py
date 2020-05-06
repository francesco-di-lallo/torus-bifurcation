# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 19:54:27 2020

@author: Francesco
"""
import matplotlib.pyplot as plt
from ode import Ode
from bmk import Bmk
import numpy as np

def para_plot(xs,ys,break_threshold):
    xs_ = []
    ys_ = []

    for index, val_pair in enumerate(zip(xs,ys)):
        if index != 0:
            if np.abs(xs[index]-xs[index-1]) > break_threshold:
                xs_.append(None)
                ys_.append(None)
            if np.abs(ys[index]-ys[index-1]) > break_threshold:
                xs_.append(None)
                ys_.append(None)
        xs_.append(val_pair[0])
        ys_.append(val_pair[1])

    return np.array(xs_), np.array(ys_)


bmk = Bmk({'ox': -0.009649961793522382, 'oy': 0.9270432434144844})
ode = Ode(bmk.ode)

P1, _ = ode.unstable_manifold(1)
P2, _ = ode.unstable_manifold(-1)

x_plt = [] 
y_plt = []
for i in range(len(P1)):
    (x,y) = ode._point_dict_getter(P1[i])
    x_plt.append(x)
    y_plt.append(y)

x_plt, y_plt = para_plot(np.array(x_plt)%1,
                         np.array(y_plt)%1, 0.5)


x_plt_lift = [] 
y_plt_lift = []
for i in range(len(P2)):
    (x,y) = ode._point_dict_getter(P2[i])
    x_plt_lift.append(x)
    y_plt_lift.append(y)

x_plt_lift, y_plt_lift = para_plot(np.array(x_plt_lift)%1,
                                   np.array(y_plt_lift)%1, 0.5)

fp = ode.get_fixed_points()
fp_list = []
for i in range(2):
    fp_list.append(ode._point_dict_getter(fp[i]))


X, Y = np.meshgrid(np.arange(0, 1, 0.1), np.arange(0, 1, 0.1))
x_shape = X.shape
U = np.zeros(x_shape)
V = np.zeros(x_shape)    

for i in range(x_shape[0]):
    for j in range(x_shape[1]):
        val = ode.evaluate(X[i,j], Y[i,j])
        U[i,j] = 0.1*val[0]
        V[i,j] = 0.1*val[1]

fig, ax = plt.subplots()
q = ax.quiver(X, Y, U, V, units='xy' , color='red')

#ax.set_aspect()
plt.plot(x_plt, y_plt)
plt.plot(x_plt_lift, y_plt_lift)
plt.scatter(*zip(*fp_list),c='k', marker = 's',zorder=3)
plt.xlim((0,1))
plt.ylim((0,1))
plt.title('ox: {}, oy: {}'.format(-0.009649961793522382, 0.9270432434144844))
plt.savefig('Necklace Point.png', dpi = 2000)
plt.show()