# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 14:19:47 2020

@author: Francesco
"""
import PyDSTool

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