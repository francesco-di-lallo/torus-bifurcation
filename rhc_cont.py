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

icdict = {'x': 0, 'y': 0} # initial conditions
pardict = {'sx': 0, 'sy': 1} # parameters

# ODE rhs
x_rhs = 'sx-cos(2*pi*y)-0.1*cos(2*pi*x)'
y_rhs = 'sy-sin(2*pi*y)-0.1*sin(2*pi*x)'

#defining ODE
vardict = {'x': x_rhs, 'y': y_rhs}

DSargs = PyDSTool.common.args()
DSargs.name = 'BMK'
DSargs.ics = icdict
DSargs.pars = pardict
DSargs.tdata = [0,25]
DSargs.varspecs = vardict
DSargs.xdomain = {'x': [0,1], 'y': [0,1]}

ode = PyDSTool.Generator.Vode_ODEsystem(DSargs)

##############################################################################

class RhcCont(PyDSTool.ContClass):
    def forward():
        pass

class Ode(object):
    
    def __init__(self, ode):
        self.ode = ode
        self.fixed_points = pp.find_fixedpoints(self.ode, n=4, eps=1e-5)
        self.in_resonance_region = True
        if len(self.fixed_points) != 2:
            print("Not in resonance region")
            self.in_resonance_region = False
    
    def _point_dict_setter(self,x,y):
        return {'x':x, 'y':y}
    
    def _point_dict_getter(self,point_dict):
        return [point_dict['x'], point_dict['y']]
    
    def get_fixed_points(self):
        return self.fixed_points
    
    def fp_eigen(self, fp_coord, fuzz_factor = 1e-5):
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
            Tolerance to zero. The default is 1e-5.

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
   
    def _initialise_unstable_manifold(self, perturbation = 1e-5):
        """
        

        Parameters
        ----------
        perturbation : float, optional
            DESCRIPTION. The default is 1e-5. Perturbation in the unstable 
            manifold direction to then make it iterable

        Returns
        -------
        saddle_coord : dict
            DESCRIPTION. Point dictionary that is now iterable via RK4 to form W-

        """
        unstable_eval_index = [i for i, mu in enumerate(self.get_saddle_point()[1]) if mu>0]                              
        unstable_evec = [self.saddle_point[2][0][unstable_eval_index],\
                         self.saddle_point[2][1][unstable_eval_index]]
        saddle_coord = self.saddle_point[0]
        saddle_coord['x'] += float(perturbation*unstable_evec[0])
        saddle_coord['y'] += float(perturbation*unstable_evec[1])
        
        return saddle_coord
            
    def _lift(self, fp):
        
        """
        Takes point dictionary and returns the lift of it in the cartesian space
        [1,2]x[0,1]

        Returns
        -------
        lift : dict

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
        norm = (x['x']-y['x'])**2+(x['y']-y['y'])**2
        return norm
    
    def evaluate(self, x,y):
        point_dict = self._point_dict_setter(x,y)
        return self.ode.Rhs(0, point_dict)
    
    def _RK4_iter(self, point_dict, h=1e-0):
        """
        Creates the next iterate of Runge-Kutta 4 method. 
        Since BKM is autonomous, it is independent of time.

        Parameters
        ----------
        point_dict : dict
            Point dictionary
        h : TYPE float, optional
            The default is 1e-5. h value in RK4 algorithm

        Returns
        -------
        list
            Next iterate of RK4. 

        """
        (x,y) = self._point_dict_getter(point_dict)
        
        (k1, l1) = self.evaluate(x,y) 
        
        (k2, l2) = self.evaluate(x+0.5*k1, y+0.5*l1)
        
        (k3, l3) = self.evaluate(x+0.5*k2, y+0.5*l2)
        
        (k4, l4) = self.evaluate(x+k3, y+l3)
        
        return self._point_dict_setter(x+(float(h)/6)*(k1+2*k2+2*k3+k4),
                y+(float(h)/6)*(l1+2*l2+2*l3+l4))

    def unstable_manifold(self):
        
        saddle = self.get_saddle_point()
        
        x1 = self._lift(saddle[0])
        p0 = self._initialise_unstable_manifold()
        P = [p0]
        pn = self._RK4_iter(p0)
        P.append(pn)
        i = 0
        while (self._dist(P[-1],x1) < self._dist(P[-2],x1)):
            P.append(self._RK4_iter(P[-1], h=0.01))
            i+=1
        
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
        return self._dist(P[0], x)
"""
Next things to be done:
    
    2) apply some grad descent in sy-param space and repeat (1), (2) 
     <finds rhc point for specific sx

    3) find a way to raise error if self.in_resonance_region = False
"""
bmk = Ode(ode)












