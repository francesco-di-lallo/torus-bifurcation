# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:34:59 2020

@author: Francesco
"""



from numpy import dot, sqrt, concatenate, zeros, cos, sin, pi
from scipy.optimize import leastsq




def f(t,x):
    """
    Rossler system written in the form of Eq. (7)
    """
    xd = zeros(len(x),'d')
    xd[0] = 0-cos(2*pi*xd[1])-0.1*cos(2*pi*xd[0])
    xd[1] = 0.92-sin(2*pi*xd[1])-0.1*sin(2*pi*xd[0])
    return xd



def integrate(t,x,func,h,w,a):
    """
    5th-order Runge-Kutta integration scheme. Input:
    t - initial time
    x - vector of initial conditions at initial time t
    h - integration time step, w - period
    a - additional parameters
    """
    k1=h*func(t,x)
    k2=h*func(t+0.5*h,x+0.5*k1)
    k3=h*func(t+0.5*h,x+(3.0*k1+k2)/16.0)
    k4=h*func(t+h,x+0.5*k3)
    k5=h*func(t+h,x+(-3.0*k2+6.0*k3+9.0*k4)/16.0)
    k6=h*func(t+h,x+(k1+4.0*k2+6.0*k3-12.0*k4+8.0*k5)/7.0)
    xp = x + (7.0*k1+32.0*k3+12.0*k4+32.0*k5+7.0*k6)/90.0
    return xp



def ef(v,x,func,dt,a,p):
    """
    Residual (error vector). Input:
    v - vector containing the quantities to be optimized
    x - vector of initial conditions
    func - function, dt - integration time step
    a - additional parameters
    p - controls length of error vector
    """
    j = int(2.0/dt)
    vv = zeros((j,len(x)),'d')
    vv[0,0:1] = v[0:1] # set initial condition
    vv[0,1] = x[1]
    T = v[1] # set period
    i = 0
    while i < j/2+p:
        t = i*dt
        vv[i+1,:] = integrate(t,vv[i,:],func,dt,T,a)
        i = i+1
    er = vv[j//2,:]-vv[0,:] # creates residual error vector
    for i in range(1,p): # of appropriate length
        er = concatenate((er,vv[j//2+i,:]-vv[i,:]))
    print('Error:', sqrt(dot(er,er)))
    return er


def main():
    a0 = zeros(3,'d') # predetermined system parameters
    a0[0] = 0.15; a0[1] = 0.2; a0[2] = 3.5
    x0 = zeros(2,'d') # initial conditions (N=3)
    x0[0] = 0.8574; x0[1] = 0.2599
    v0 = zeros(2,'d') # quantities for optimization
    v0[0:1] = x0[0:1]
    v0[1] = 5.92030065 # initial guess for period
    p = 2 # length of residual is pN
    h = 1.0/1024.0 # integration time step
    # # LM optimization
    v, succ = leastsq(ef,v0,args=(x0,f,h,a0,p),ftol=1e-12,maxfev=200)
    err = ef(v,x0,f,h,a0,p) # error estimation
    es = sqrt(dot(err,err))
    #
    u0 = (v[0],v[1],x0[2],v[1],es/1e-13)
    print ('IC: %20.16f %20.16f %20.16f Period: %20.16f Magnitude %6.2f' % u0)

main()