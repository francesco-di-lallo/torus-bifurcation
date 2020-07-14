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
        """
        Instantiate RHC continuation class

        """
        print("Computing RHC Curve")
        
    
    def learning(self, ox, oy_seed, direction):
        """
        Learns the zero of PE for a given ox parameter for a given branch

        Parameters
        ----------
        ox : float
            DESCRIPTION. Fixed ox parameter 
        oy_seed : float
            DESCRIPTION. Seed for NR method
        direction : +/- 1
            DESCRIPTION. +1 for positive branch of unstable manifold. 
            -1 for negative

        Returns
        -------
        estimate : float
            DESCRIPTION. Estimate of the RHC point given ox. To within 1e-5.

        """
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
        """
        Computes an iterate of NR

        Parameters
        ----------
        ox : float
            DESCRIPTION. 
        oy : float
            DESCRIPTION.
        direction : +/- 1
            DESCRIPTION. 1 for positive branch of unstable manifold. 
            -1 for negative
        step : TYPE, optional
            DESCRIPTION. The default is 1e-5.

        Returns
        -------
        PE/PE' : float
            DESCRIPTION. Error iterate

        """
        bmk = Bmk({'ox':ox, 'oy':oy})
        ode = Ode(bmk.ode)
        
        PE = ode.homoclinic_dist(direction, 'rotational')
        
        
        bmk_incr = Bmk({'ox':ox, 'oy':oy+step})
        ode_incr = Ode(bmk_incr.ode)
        
        PE_incr = ode_incr.homoclinic_dist(direction, 'rotational')
        return PE/((PE_incr-PE)/step)
    
    def seed(self, direction):
        """
        Computes the seed of a given direction with seed (0,1)

        Parameters
        ----------
        direction : +/- 1
            DESCRIPTION. 1 for positive branch of unstable manifold. 
            -1 for negative

        Returns
        -------
        oy seed : oy value for ox = 0 on the given branch of RHC.
            DESCRIPTION.

        """
        return self.learning(0, 1, direction)
    
    
    def propagate(self, direction):
        """
        Propagates the initial seed of RHC on a given branch

        Parameters
        ----------
        direction : +/- 1
            DESCRIPTION. 1 for positive branch of unstable manifold. 
            -1 for negative

        Returns
        -------
        rhc_ox : list
            DESCRIPTION. ox values of RHC curve
        rhc_oy : list
            DESCRIPTION. oy values of RHC curve

        """
        
        if direction == 1:
            rhc_oy = pickle.load(open('rhc_oy_heun.p', 'rb'))
            rhc_ox = pickle.load(open('rhc_ox_heun.p', 'rb'))
        else:
            rhc_oy = pickle.load(open('rhc_oy_lift_heun.p', 'rb'))
            rhc_ox = pickle.load(open('rhc_ox_lift_heun.p', 'rb'))
        """
        print('Initialising Seed')
        rhc_ox = [0]
        rhc_oy = [self.seed(direction)]
        """
        ox_step = 0.0001
        print("Propagating RHC Forwards")
        for i in range(0):
            print("{0}% Iterating {1} Forwards".format(100*i/20, direction))
            ox_seed = rhc_ox[-1]+ox_step
            oy_seed = rhc_oy[-1]
            try:
                rhc_oy.append(self.learning(ox_seed, oy_seed, direction))
                rhc_ox.append(ox_seed)
            except:
                print("done goofed")
            
            if direction == 1:
                pickle.dump(rhc_ox, open('rhc_ox_heun.p', 'wb'))
                pickle.dump(rhc_oy, open('rhc_oy_heun.p', 'wb'))
            else:
                pickle.dump(rhc_ox, open('rhc_ox_lift_heun.p', 'wb'))
                pickle.dump(rhc_oy, open('rhc_oy_lift_heun.p', 'wb')) 
              
        print("Propagating RHC Backwards")        
        for i in range(50):
            print("{0}% Iterating {1} Backwards".format(100*i/10., direction))
            ox_seed = rhc_ox[0]-ox_step
            oy_seed = rhc_oy[0]
            try:
                rhc_oy.insert(0, self.learning(ox_seed, oy_seed, direction))
                rhc_ox.insert(0, ox_seed)
            except:
                print("done goofed")
            
            if direction == 1:
                pickle.dump(rhc_ox, open('rhc_ox_heun.p', 'wb'))
                pickle.dump(rhc_oy, open('rhc_oy_heun.p', 'wb'))
            else:
                pickle.dump(rhc_ox, open('rhc_ox_lift_heun.p', 'wb'))
                pickle.dump(rhc_oy, open('rhc_oy_lift_heun.p', 'wb')) 
            
        return rhc_ox, rhc_oy


RHC = RhcCont()
#RHC.propagate(1)
RHC.propagate(-1)