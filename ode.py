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
class Ode(object):
    
    def __init__(self, ode):
        """
        Saves ODE as self.ode
        Sets fixed points as self.fixed_points
        Verifies if the parameters are in the resonance region

        Parameters
        ----------
        ode : PyDSTool.Generator.Vode_ODEsystem.Vode_ODEsystem
            DESCRIPTION. BMK type ODE
        """
        
        
        self.ode = ode
        self.fixed_points = pp.find_fixedpoints(self.ode, n=10, eps=1e-5)
        self.in_resonance_region = True
        if len(self.fixed_points) != 2:
            print("Not in resonance region")
            self.in_resonance_region = False
            
    
    def _point_dict_setter(self,x,y):
        """
        Takes two coordinates and sets it to a point dictionary

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
        """
        Gets fixed points as point dictionary

        """
        return self.fixed_points
    
    def fp_eigen(self, fp_coord, fuzz_factor = 1e-12):
        """
        Returns information about the fixed point:
            - Stability type
            - Eigenvalues and corresponding eigenvectors

        Parameters
        ----------
        fp_coord : dict
           DESCRIPTION.  Dictionary {'x': x-coordinate of fixed point,
                        'y': y-coordinate of fixed point.
        fuzz_factor : float, optional
            DESCRIPTION. Tolerance to zero. The default is 1e-12.

        Returns
        -------
        fp_type : str
           DESCRIPTION.  Description of type of fixed point:
                - saddle point
                - hyperbolic attractor
                - hyperbolic repeller
                - non-hyperbolic
        eigen_vectors : array
           DESCRIPTION.  First entry contains eigenvalues
            Second entry contains eigenvectors (v1|v2)

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
            fp_type = "centre"
    
        if fp_type != "non-hyperbolic":
            return (fp_type, eigen_vectors)
        else:
            return (fp_type, None)
        
    def jacobian(self, fp):
        """
        Parameters
        ----------
        fp : dict
           DESCRIPTION.  Dictionary of a point

        Returns
        -------
        j : 2x2 array
           DESCRIPTION.  Computed jacobian at the given point.
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
            DESCRIPTION. First entry is a dictionary of the coordinates of the saddle point
            Second entry is a 2x1 array of the eigenvalues [lambda_1, lambda_2]
            Third entry is a 2x2 array of the eigenvectors [v1|v2]

        """

        if self.in_resonance_region:
            for fp in self.fixed_points:
                if self.fp_eigen(fp)[0] == 'saddle point':
                        self.saddle_point = (fp, self.fp_eigen(fp)[1][0],\
                                             self.fp_eigen(fp)[1][1])
                        return self.saddle_point
    def get_centre(self):
        if self.in_resonance_region:
            for fp in self.fixed_points:
                if self.fp_eigen(fp)[0] in ['hyperbolic attractor',\
                                            'hyperbolic repeller', 'centre']:
                        self.centre_point = (fp, self.fp_eigen(fp)[1][0],\
                                             self.fp_eigen(fp)[1][1])
                        return self.centre_point
   
    def _initialise_unstable_manifold(self,direction, perturbation = 1e-7):
        """
        Returns a point dictionary that is perturbed along the unstable eigenvector.
        This is a linear approximation of the unstable manifold and 
        is then used as the ic of the trajectory to compute W^-

        Parameters
        ----------
        direction : -1 or +1
            DESCRIPTION. +1 corresponds to unstable branch of saddle. -1 corresponds to 
            stable branch.
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
            
        saddle_coord = self.get_saddle_point()[0].copy()
        
        if direction == -1:
            saddle_coord['x'] += 1
            
        if np.sign(unstable_evec[0]) == direction:
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
            DESCRIPTION. L2 norm of x, y.

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

    def _euler_iter(self, point_dict, h=1e-6):
        """
        Computes the next Euler-Heun iterate

        Parameters
        ----------
        point_dict : Point dictionary 
            DESCRIPTION. x_n
        h : float, optional
            DESCRIPTION. The default is 1e-6.

        Returns
        -------
        x_n+1 : Point dictionary

        """
        
        (x,y) = self._point_dict_getter(point_dict)
        
        (x_dot, y_dot) = self.evaluate(x,y)
        
        (x_tilde, y_tilde) = (x+float(h)*x_dot, y+float(h)*y_dot)
        
        (x_tilde_dot, y_tilde_dot) = self.evaluate(x_tilde,y_tilde)
        
        (x_next, y_next) = (x+float(h/2)*(x_dot+x_tilde_dot),
                            y+float(h/2)*(y_dot+y_tilde_dot))
        
        return self._point_dict_setter(x_next, y_next)
        
    def unstable_manifold(self, direction):
        """
        Computes the branch of the unstable manifold in the specified direction
        Computed until the distance to the non-starting saddle increases.

        Parameters
        ----------
        direction : +/- 1
            DESCRIPTION. +1 corresponds to unstable branch of saddle. -1 corresponds to 
            stable branch.

        Returns
        -------
        P :  List of Point dictionaries
            DESCRIPTION. Point dictionaries forming the manifold
        x1 : TYPE
            DESCRIPTION. The saddle which the manifold is approaching

        """
        
        saddle = self.get_saddle_point()
        if direction == 1:
            x1 = self._lift(saddle[0])
        else:
            x1 = saddle[0]
        p0 = self._initialise_unstable_manifold(direction, 1e-3) ##
        P = [p0]
        p1 = self._euler_iter(p0)
        P.append(p1)
        i = 0
        
        while (self._dist(P[-1],x1) < self._dist(P[-2],x1)) or i < 500:
            P.append(self._euler_iter(P[-1], h=1e-3))
            i += 1
        
        return P, x1
    
    def extended_unstable_manifold(self, direction):
        P, x1 = self.unstable_manifold(direction)
        
        x0 = P[0]
        i = 0
        while (self._dist(P[-1],x0) < self._dist(P[-2],x0)) or i < 500:
            P.append(self._euler_iter(P[-1], h=1e-3))
            i += 1
        return P, x0
        
        
    def plot(self, dict_list):
        """
        Implementing matplotlib.pyplot.plot but passing an array of
        point dictionaries instead of an array of floats.

        Parameters
        ----------
        dict_list : Point dictionary

        """
        x_plt = [] 
        y_plt = []
        for i in range(len(dict_list)):
            (x,y) = self._point_dict_getter(dict_list[i])
            x_plt.append(x)
            y_plt.append(y)
        
        plt.plot(x_plt,y_plt)
        plt.show()
        
    def homoclinic_dist(self, direction, t):
        """
        Computes the Pontryagin energy for a given direction

        Parameters
        ----------
        direction : +/- 1 
            DESCRIPTION. +1 corresponds to unstable branch of saddle. -1 corresponds to 
            stable branch.

        Returns
        -------
        Pontryagin Energy : float
            DESCRIPTION. PE+/-

        """
        if t == 'rotational':
            P, x = self.unstable_manifold(direction)
            return self._dist(P[-1], x)
        elif t == 'contractible':
            P, x = self.extended_unstable_manifold(direction)
            return self._dist(P[-1], x)
    
    def pontryagin_energy(self, t):
        """
        Computes the tuple (PE+, PE-)
        
        Parameters
        ----------
        t : 'rotational' or 'contractible'

        Returns
        -------
        PE : tuple
            DESCRIPTION. (PE+, PE-)

        """
        return (self.homoclinic_dist(1, t), 
                self.homoclinic_dist(-1, t))
    
        