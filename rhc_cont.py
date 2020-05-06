# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 14:20:59 2020

@author: Francesco
"""
from ode import Ode
from bmk import Bmk
import numpy as np
import pickle

class RhcCont(object):
    def __init__(self):
        print("Computing RHC Curve")
        
    
    def learning(self, ox, oy_seed, direction):
        oy = oy_seed 
        learning_rate = 1e-1
        estimates = [oy_seed]
        oy -= learning_rate*self.newton_raphson(ox, oy, direction)
        estimates.append(oy)
        i = 0
        while np.abs(estimates[-1]-estimates[-2])>1e-5 and i<100:
            print(oy)
            oy -= learning_rate*self.newton_raphson(ox, oy, direction)
            estimates.append(oy)
            i+=1
        return estimates[-1]
    
    def newton_raphson(self, ox, oy, direction, step = 1e-5):
        bmk = Bmk({'ox':ox, 'oy':oy})
        ode = Ode(bmk.ode)
        
        PE = ode.homoclinic_dist(direction)
        
        
        bmk_incr = Bmk({'ox':ox, 'oy':oy+step})
        ode_incr = Ode(bmk_incr.ode)
        
        PE_incr = ode_incr.homoclinic_dist(direction)
        return PE/((PE_incr-PE)/step)
    
    def seed(self, direction):
        return self.learning(0, 1, direction)
    
    
    def propagate(self, direction):
        """
        if direction == 1:
            rhc_oy = pickle.load(open('rhc_oy_05.p', 'rb'))
            rhc_ox = pickle.load(open('rhc_ox_05.p', 'rb'))
        else:
            rhc_oy = pickle.load(open('rhc_oy_lift_05.p', 'rb'))
            rhc_ox = pickle.load(open('rhc_ox_lift_05.p', 'rb'))
        """
        print('Initialising Seed')
        rhc_ox = [0]
        rhc_oy = [self.seed(direction)]
        
        ox_step = 0.01
        print("Propagating RHC Forwards")
        for i in range(50):
            print("{0}% Iterating {1} Forwards".format(i/4., direction))
            ox_seed = rhc_ox[-1]+ox_step
            oy_seed = rhc_oy[-1]
            try:
                rhc_oy.append(self.learning(ox_seed, oy_seed, direction))
                rhc_ox.append(ox_seed)
            except:
                print("done goofed")
            """
            if direction == 1:
                pickle.dump(rhc_ox, open('rhc_ox_05.p', 'wb'))
                pickle.dump(rhc_oy, open('rhc_oy_05.p', 'wb'))
            else:
                pickle.dump(rhc_ox, open('rhc_ox_lift_05.p', 'wb'))
                pickle.dump(rhc_oy, open('rhc_oy_lift_05.p', 'wb')) 
            """    
        print("Propagating RHC Backwards")        
        for i in range(50):
            print("{0}% Iterating {1} Backwards".format(i/4., direction))
            ox_seed = rhc_ox[0]-ox_step
            oy_seed = rhc_oy[0]
            try:
                rhc_oy.insert(0, self.learning(ox_seed, oy_seed, direction))
                rhc_ox.insert(0, ox_seed)
            except:
                print("done goofed")
            """
            if direction == 1:
                pickle.dump(rhc_ox, open('rhc_ox_05.p', 'wb'))
                pickle.dump(rhc_oy, open('rhc_oy_05.p', 'wb'))
            else:
                pickle.dump(rhc_ox, open('rhc_ox_lift_05.p', 'wb'))
                pickle.dump(rhc_oy, open('rhc_oy_lift_05.p', 'wb')) 
            """
        return rhc_ox, rhc_oy
