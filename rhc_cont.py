# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 18:26:28 2020

@author: Francesco
"""


import PyDSTool
from numpy import (sin, cos, pi)
import numpy as np
from PyDSTool.Toolbox import phaseplane as pp
import matplotlib.pyplot as plt
from bmk import Bmk

##############################################################################

class RhcCont(object):
    
    
    def learning(self, ox, oy_seed):
        oy = oy_seed 
        learning_rate = 1e-4
        estimates = [oy_seed]
        oy -= learning_rate*self.grad_PE(ox, oy)
        estimates.append(oy)
        while np.abs(estimates[-1]-estimates[-2])>1e-4:
            oy -= learning_rate*self.grad_PE(ox, oy)
            estimates.append(oy)
        return estimates[-1]
    
    def grad_PE(self, ox, oy, step = 0.01):
        bmk = Bmk({'ox':ox, 'oy':oy})
        ode = Ode(bmk.ode)
        
        PE = ode.homoclinic_dist()
        
        
        bmk_incr = Bmk({'ox':ox, 'oy':oy+step})
        ode_incr = Ode(bmk_incr.ode)
        
        PE_incr = ode_incr.homoclinic_dist()
        return (PE_incr-PE)/step
    def seed(self):
        return self.learning(0, 0.91)

    def propagate_forwards(self):
        rhc_oy = [self.seed()]
        rhc_ox = [0]
        self.seed = rhc_oy[0]
        ox_step = 0.05
        for i in range(10):
            print(i)
            ox_seed = rhc_ox[-1]+ox_step
            oy_seed = rhc_oy[-1]
            try:
                rhc_ox.append(ox_seed)
                rhc_oy.append(self.learning(ox_seed, oy_seed))
            except:
                print("done goofed")
        return rhc_ox, rhc_oy

class Ode(object):
    
    def __init__(self, ode):
        """
        Saves ODE as self.ode
        Sets fixed points as self fixed_points
        Verifies if the parameters are in the resonance region

        Parameters
        ----------
        ode : PyDSTool.Generator.Vode_ODEsystem.Vode_ODEsystem
            DESCRIPTION. BMK type ODE

        Returns
        -------
        None.

        """
        self.ode = ode
        self.fixed_points = pp.find_fixedpoints(self.ode, n=10, eps=1e-5)
        self.in_resonance_region = True
        if len(self.fixed_points) != 2:
            print("Not in resonance region")
            self.in_resonance_region = False
    
    def _point_dict_setter(self,x,y):
        """
        Takes a cartesian point and sets it to a point dictionary

        Parameters
        ----------
        x : float
            DESCRIPTION. x coordinate
        y : float
            DESCRIPTION. y coordinate

        Returns
        -------
        dict

        """
        return {'x':x, 'y':y}
    
    def _point_dict_getter(self,point_dict):
        """
        Converts point dictionary to a list

        Parameters
        ----------
        point_dict : dict

        Returns
        -------
        list

        """
        return [point_dict['x'], point_dict['y']]
    
    def get_fixed_points(self):
        return self.fixed_points
    
    def fp_eigen(self, fp_coord, fuzz_factor = 1e-12):
        """
        Returns information about the fixed point:
            - Stability type
            - Eigenvalues and corresponding eigenvectors

        Parameters
        ----------
        fp_coord : dict
            Dictionary {'x': x-coordinate of fixed point,
                        'y': y-coordinate of fixed point.
        fuzz_factor : float, optional
            Tolerance to zero. The default is 1e-12.

        Returns
        -------
        fp_type : str
            Description of type of fixed point:
                - saddle point
                - hyperbolic attractor
                - hyperbolic repeller
                - non-hyperbolic
        eigen_vectors : array
            First array contains eigenvalues
            Second array contains eigenvectors (v1|v2)

        """
        jac = self.jacobian(fp_coord)
        eigen_values = np.linalg.eigvals(jac)
        eigen_vectors = np.linalg.eig(jac)
        fp_type = 'unknown'
    
        # stability
        det_jac = np.linalg.det(jac)
        tr_jac = np.trace(jac)
    
        if abs(det_jac) < fuzz_factor:
            fp_type = "non-hyperbolic"
        elif np.sometrue(np.abs(np.real(eigen_values)) > fuzz_factor):
            if det_jac < 0:
                fp_type = "saddle point"
            if det_jac > 0:                
                    if tr_jac < 0:
                        fp_type = "hyperbolic attractor"
                    if tr_jac > 0:
                        fp_type = "hyperbolic repeller"
        else:
            fp_type = "non-hyperbolic"
    
        if fp_type != "non-hyperbolic":
            return (fp_type, eigen_vectors)
        else:
            return (fp_type, None)
        
    def jacobian(self, fp):
        """
        Parameters
        ----------
        fp : dict
            Dictionary of a point

        Returns
        -------
        j : 2x2 array
            Computed jacobian at the given point.
            Note that the jacobian is independent of the parameter
            so all BMK ODEs have the same jacobian

        """
        (x,y) = self._point_dict_getter(fp)
        j = [[2*pi*0.1*sin(2*pi*x), 2*pi*sin(2*pi*y)],
             [-2*pi*0.1*cos(2*pi*x), -2*pi*cos(2*pi*y)]]
        return j

    def get_saddle_point(self):
        """
        Getter/setter method that retrieves information on the saddle fixed point
        provided that the ODE is in the resonance region

        Returns 
        -------
        tuple
            First entry is a dictionary of the coordinates of the saddle point
            Second entry is a 2x1 array of the eigenvalues [lambda_1, lambda_2]
            Third entry is a 2x2 array of the eigenvectors [v1|v2]

        """

        if self.in_resonance_region:
            for fp in self.fixed_points:
                if self.fp_eigen(fp)[0] == 'saddle point':
                        self.saddle_point = (fp, self.fp_eigen(fp)[1][0],\
                                             self.fp_eigen(fp)[1][1])
                        return self.saddle_point
   
    def _initialise_unstable_manifold(self, perturbation = 1e-7):
        """
        Sets a point dictionary that is perturbed along the unstable eigenvector.
        This is a linear approximation of the unstable manifold and 
        is then used as the ic of the trajectory to compute W^-

        Parameters
        ----------
        perturbation : float, optional
            DESCRIPTION. The default is 1e-5. Perturbation in the unstable 
            manifold direction to then make it iterable.

        Returns
        -------
        saddle_coord : dict
            DESCRIPTION. Point dictionary near saddle point.

        """
        unstable_eval_index = [i for i, mu in enumerate(self.get_saddle_point()[1]) if mu>0]                              
        unstable_evec = [self.saddle_point[2][0][unstable_eval_index],\
                         self.saddle_point[2][1][unstable_eval_index]]
        saddle_coord = self.saddle_point[0]
        if unstable_evec[0]>0:
            saddle_coord['x'] += float(perturbation*unstable_evec[0])
            saddle_coord['y'] += float(perturbation*unstable_evec[1])
        else:
            saddle_coord['x'] -= float(perturbation*unstable_evec[0])
            saddle_coord['y'] -= float(perturbation*unstable_evec[1])
        
        return saddle_coord
            
    def _lift(self, fp):
        """
        Takes point dictionary and returns the lift of it in the cartesian space
        [1,2]x[0,1]

        Parameters
        ----------
        fp : dict
            DESCRIPTION. Point dict 

        Returns
        -------
        lift : dict
            DESCRIPTION. Point dict

        """

        lift = {'x': fp['x']+1,
                'y': fp['y']}
        return lift

    def _dist(self,x,y):
        """
        dist is the L2 distance between x and y

        Parameters
        ----------
        x : point dict
            
        y : point dict
            

        Returns
        -------
        norm : float
            L2 norm of x, y.

        """
        norm = np.abs(x['x']-y['x'])**2+np.abs(x['y']-y['y'])**2
        return norm
    
    def evaluate(self, x,y):
        """
        Evaluates the ODE as a vector field.

        Parameters
        ----------
        x : float

        y : float

        Returns
        -------
        array

        """
        point_dict = self._point_dict_setter(x,y)
        return self.ode.Rhs(0, point_dict)

    def _euler_iter(self, point_dict, h=1e-6): ##
        
        (x,y) = self._point_dict_getter(point_dict)
        
        (x_dot, y_dot) = self.evaluate(x,y)
        
        return self._point_dict_setter(x+float(h)*x_dot, y+float(h)*y_dot)
        
    def unstable_manifold(self):
        
        saddle = self.get_saddle_point()
        
        x1 = self._lift(saddle[0])
        p0 = self._initialise_unstable_manifold(1e-3) ##
        P = [p0]
        p1 = self._euler_iter(p0)
        P.append(p1)
        i = 0
        while (self._dist(P[-1],x1) < self._dist(P[-2],x1)) or i < 500:
            P.append(self._euler_iter(P[-1], h=1e-3))
            i += 1
        
        return P, x1
    def plot(self, dict_list):
        x_plt = [] 
        y_plt = []
        for i in range(len(dict_list)):
            (x,y) = self._point_dict_getter(dict_list[i])
            x_plt.append(x)
            y_plt.append(y)
        
        plt.plot(x_plt,y_plt)
        plt.show()
    def homoclinic_dist(self):
        P, x = self.unstable_manifold()
        return self._dist(P[-1], x)
    


dist = []
oy_list = []

for oy in np.arange(0.9,1, 0.001):
    bmk = Bmk({'ox':0.5, 'oy':oy})
    ode = Ode(bmk.ode)
    print(oy)
    if ode.in_resonance_region:
        dist.append(ode.homoclinic_dist())
        oy_list.append(oy)

plt.plot(oy_list, dist)
plt.title('Pontryagin Energy for ox = 0.5')
plt.xlabel('oy')
plt.ylabel('Pontryagin Energy')
plt.xlim=(0.9,1)
plt.ylim=(0,1)

plt.show()

"""
Next things to be done:

    3) To compute other RHC: run unstable manifold but iterate from the lift 
        of the saddle in the negative perturbed direction
"""