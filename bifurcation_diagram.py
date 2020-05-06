# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:59:47 2020

@author: Francesco
"""

import PyDSTool
#from numpy import (sin, cos, pi)
#from PyDSTool.Toolbox import phaseplane as pp
#import numpy as np
import matplotlib.pyplot as plt
from rhc_cont import Ode
import pickle
from bmk import Bmk
#from ode import Ode
from rhc_cont import RhcCont


bmk = Bmk({'ox': 0, 'oy': 0})
ode = bmk.get_ode()

#Setting Continuation class
PC = PyDSTool.ContClass(ode)
##########################################################
#
#
#
##########################################################
def inner_outer_init(PC):
    print("Initialising boundary")
    PCargs = PyDSTool.args(name='EQ1', type='EP-C')
    PCargs.freepars     = ['ox']
    PCargs.MaxNumPoints = 200       
    PCargs.MaxStepSize  = 0.001
    PCargs.LocBifPoints = 'LP'
    
    PC.newCurve(PCargs)
    PC['EQ1'].forward()
    return PC
##########################################################
#
#
#
##########################################################
def inner_sne(PC):
    print("Computing Inner SNE")
    PCargs = PyDSTool.args(name='SNE1', type='LP-C')
    PCargs.initpoint = 'EQ1:LP1'
    PCargs.freepars  = ['ox', 'oy']
    PCargs.MaxStepSize = 0.1
    PCargs.LocBifPoints = ['BT']
    PCargs.MaxNumPoints = 100
    PC.newCurve(PCargs)
    PC['SNE1'].forward()
    return PC
##########################################################
#
#
#
##########################################################
def outer_sne(PC):
    print("Computing Outer SNE")
    PCargs = PyDSTool.args(name='SNE2', type='LP-C')
    PCargs.initpoint = 'EQ1:LP2'
    PCargs.freepars  = ['ox', 'oy']
    PCargs.MaxStepSize = 0.1
    PCargs.LocBifPoints = ['BT']
    PCargs.MaxNumPoints = 100
    PC.newCurve(PCargs)
    PC['SNE2'].forward()
    return PC
##########################################################
#
#
#
##########################################################
def trace_zero_upper(PC):
    print("Computing tr0 loop upper")
    PCargs = PyDSTool.args(name='TR01', type='H-C2')
    PCargs.initpoint = 'SNE1:BT1'
    PCargs.freepars  = ['ox', 'oy']
    PCargs.MaxStepSize = 0.1
    PCargs.LocBifPoints = ['BT']
    PCargs.MaxNumPoints = 180
    PC.newCurve(PCargs)
    PC['TR01'].forward()
    return PC
##########################################################
#
#
#
##########################################################
def trace_zero_lower(PC):
    print("Computing tr0 loop lower")
    PCargs = PyDSTool.args(name='TR02', type='H-C2')
    PCargs.initpoint = 'SNE1:BT2'
    PCargs.freepars  = ['ox', 'oy']
    PCargs.MaxStepSize = 0.1
    PCargs.LocBifPoints = ['BT']
    PCargs.MaxNumPoints = 180
    PC.newCurve(PCargs)
    PC['TR02'].forward()
    return PC
##########################################################
#
#
#
##########################################################
def bifurcation_diagram(PC):
    
    inner_outer_init(PC)
    
    inner_sne(PC)
    outer_sne(PC)
    trace_zero_lower(PC)
    trace_zero_upper(PC)
    
    PC['SNE1'].display(coords = ['ox','oy'], linewidth = 0.2)
    
    PC['SNE2'].display(coords = ['ox','oy'], linewidth = 0.2)
    
    PC['TR01'].display(coords = ['ox','oy'], linewidth = 0.2)
    
    PC['TR02'].display(coords = ['ox','oy'], linewidth = 0.2)
    
    PC.plot.togglePoints(visible='off')
    PC.plot.toggleLabels(visible='off')
    
    RHC = RhcCont()
    
    rhc_oy = pickle.load(open('rhc_oy.p', 'rb'))
    rhc_ox = pickle.load(open('rhc_ox.p', 'rb'))
    rhc_oy_lift = pickle.load(open('rhc_oy_lift.p', 'rb'))
    rhc_ox_lift = pickle.load(open('rhc_ox_lift.p', 'rb'))
    """
    rhc_ox, rhc_oy = RHC.propagate(1)
    rhc_ox_lift, rhc_oy_lift = RHC.propagate(-1)
    """
    plt.plot(rhc_ox, rhc_oy, linewidth = 0.2)
    plt.plot(rhc_ox_lift, rhc_oy_lift, linewidth = 0.2)
    
    bottom_rhc_ox = [-1*ox for ox in rhc_ox]
    bottom_rhc_oy = [-1*oy for oy in rhc_oy]
    bottom_rhc_ox_lift = [-1*ox for ox in rhc_ox_lift]
    bottom_rhc_oy_lift = [-1*oy for oy in rhc_oy_lift]
    
    plt.plot(bottom_rhc_ox, bottom_rhc_oy, linewidth = 0.2)
    plt.plot(bottom_rhc_ox_lift, bottom_rhc_oy_lift, linewidth = 0.2)
    
    plt.scatter(rhc_ox[0], rhc_oy[0], s = 0.02, c='k', marker = 's',zorder=3)
    plt.scatter(rhc_ox[-1], rhc_oy[-1], s = 0.02, c='k', marker = 's',zorder=3)
    plt.scatter(rhc_ox_lift[0], rhc_oy_lift[0], s = 0.02, c='k', marker = 's',zorder=3)
    plt.scatter(rhc_ox_lift[-1], rhc_oy_lift[-1], s = 0.02, c='k', marker = 's',zorder=3)
    plt.scatter([-0.009649961793522382], [0.9270432434144844], s = 0.01, c='k', marker = 's',zorder=3)
    
    plt.scatter(-1*rhc_ox[0], -1*rhc_oy[0], s = 0.02, c='k', marker = 's',zorder=3)
    plt.scatter(-1*rhc_ox[-1], -1*rhc_oy[-1], s = 0.02, c='k', marker = 's',zorder=3)
    plt.scatter(-1*rhc_ox_lift[0], -1*rhc_oy_lift[0], s = 0.02, c='k', marker = 's',zorder=3)
    plt.scatter(-1*rhc_ox_lift[-1], -1*rhc_oy_lift[-1], s = 0.02, c='k', marker = 's',zorder=3)
    plt.scatter([0.009649961793522382], [-0.9270432434144844], s = 0.02, c='k', marker = 's',zorder=3)
    
    plt.title('Bifurcation Diagram')    
    plt.savefig('Bifurcation Diagram.png', dpi = 2000)
    
    plt.xlim(-0.15, 0.15)
    plt.ylim(0.87,1.12)
    plt.savefig('Bifurcation Diagram Zoom.png', dpi = 2000)
    plt.show()
    
if __name__ == '__main__':
    bifurcation_diagram(PC)