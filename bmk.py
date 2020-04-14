# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:59:47 2020

@author: Francesco
"""

import PyDSTool
from numpy import (sin, cos, pi)
from PyDSTool.Toolbox import phaseplane as pp


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

class Bmk(object):
    def __init__(self, pardict):
        icdict = {'x': 0, 'y': 0} # initial conditions
        
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
        DSargs.xdomain = {'x': [0, 1], 'y': [0, 1]}
        DSargs.varspecs = vardict
        
        self.ode = PyDSTool.Generator.Vode_ODEsystem(DSargs)
        
    def get_ode(self):
        return self.ode
        

class RHCContClass(PyDSTool.ContClass):
        
    def forward(self):
        pass

#Setting Continuation class
PC = PyDSTool.ContClass(ode)
##########################################################
#
#
#
##########################################################
def inner_outer_init(PC):
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
    PCargs = PyDSTool.args(name='SNE2', type='LP-C')
    PCargs.initpoint = 'EQ1:LP2'
    PCargs.freepars  = ['ox', 'oy']
    PCargs.MaxStepSize = 0.1
    PCargs.LocBifPoints = ['BT']
    PCargs.MaxNumPoints = 80
    PC.newCurve(PCargs)
    PC['SNE2'].forward()
    return PC
##########################################################
#
#
#
##########################################################
def trace_zero_upper(PC):
    PCargs = PyDSTool.args(name='NS1', type='H-C2')
    PCargs.initpoint = 'SNE1:BT1'
    PCargs.freepars  = ['ox', 'oy']
    PCargs.MaxStepSize = 0.2
    PCargs.LocBifPoints = ['BT']
    PCargs.MaxNumPoints = 100
    PC.newCurve(PCargs)
    PC['NS1'].forward()
    return PC
##########################################################
#
#
#
##########################################################
def trace_zero_lower(PC):
    PCargs = PyDSTool.args(name='NS2', type='H-C2')
    PCargs.initpoint = 'SNE1:BT2'
    PCargs.freepars  = ['ox', 'oy']
    PCargs.MaxStepSize = 0.2
    PCargs.LocBifPoints = ['BT']
    PCargs.MaxNumPoints = 100
    PC.newCurve(PCargs)
    PC['NS2'].forward()
    return PC
##########################################################
#
#
#
##########################################################
"""
inner_outer_init(PC)
inner_sne(PC)
outer_sne(PC)
trace_zero_lower(PC)
trace_zero_upper(PC)
"""
#PC.display(['ox','oy'])