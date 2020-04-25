# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:59:47 2020

@author: Francesco
"""

import PyDSTool
from numpy import (sin, cos, pi)
from PyDSTool.Toolbox import phaseplane as pp
import numpy as np
import matplotlib.pyplot as plt
from rhc_cont import Ode
import pickle
from bmk import Bmk
from ode import Ode
from rhc_cont import RhcCont


icdict = {'x': 0, 'y': 0} # initial conditions
pardict = {'ox': -2, 'oy': 0} # parameters

# ODE rhs
x_rhs = 'ox-cos(2*pi*y)-0.1*cos(2*pi*x)'
y_rhs = 'oy-sin(2*pi*y)-0.1*sin(2*pi*x)'

#defining ODE
vardict = {'x': x_rhs, 'y': y_rhs}

DSargs = PyDSTool.common.args()
DSargs.name = 'BMK'
DSargs.ics = icdict
DSargs.pars = pardict
DSargs.tdata = [0,25]
DSargs.varspecs = vardict

ode = PyDSTool.Generator.Vode_ODEsystem(DSargs)

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
    PCargs.MaxNumPoints = 120        
    PCargs.MaxStepSize  = 0.01
    PCargs.MinStepSize  = 0.001
    PCargs.StepSize     = 0.005
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
    PCargs.MaxNumPoints = 80
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
    PCargs.MaxStepSize = 0.2
    PCargs.LocBifPoints = ['BT']
    PCargs.MaxNumPoints = 110
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
    PCargs.MaxStepSize = 0.2
    PCargs.LocBifPoints = ['BT']
    PCargs.MaxNumPoints = 110
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
    PC.plot.setLabels('BT', bytype='BT')
    
    RHC = RhcCont()
    rhc_ox, rhc_oy = RHC.propagate(1)
    rhc_ox_lift, rhc_oy_lift = RHC.propagate(-1)
    plt.plot(rhc_ox, rhc_oy, linewidth = 0.2)
    plt.plot(rhc_ox_lift, rhc_oy_lift, linewidth = 0.2)
    
    plt.title('Bifurcation Diagram')
    plt.savefig('Bifurcation Diagram.png', dpi = 2000)
    plt.show()