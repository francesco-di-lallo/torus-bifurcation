# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 18:26:22 2020

@author: Francesco
"""
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

eps = 1e-1

bound = np.arctan(eps-np.sqrt(eps**2+1))/np.pi
t_tr0_centre = np.linspace(-0.15, 0.1,10000)#bound,-bound,1000)

print('Computing P and Q')
t, x, y = sp.symbols('t x y')

S = sp.sin(2*np.pi*t)
C = sp.cos(2*np.pi*t)

Gamma = np.sqrt(1-(eps**2) * (S**2))

f = np.sqrt(eps*(Gamma*C-eps*S**2))

omega = 2*np.pi*f

a = eps*S*(Gamma+eps*C)
b = 1
c = f*Gamma
d = -f*eps*S

Delta = a*d-b*c

P = (2*np.pi**2/(Delta**2))*((a*eps*C+b*eps*S)*(d*x-b*y)**2 + \
                             (a*eps*S+b*Gamma)*(-c*x+a*y)**2) +\
    (4*np.pi**3/(3*Delta**3))*((-a*eps*S+b*eps*C)*(d*x-b*y)**3 + \
                             (b*eps*S-a*Gamma)*(-c*x+a*y)**3)
        
Q = (2*np.pi**2/(Delta**2))*((c*eps*C+d*eps*S)*(d*x-b*y)**2 + \
                             (c*eps*S+d*Gamma)*(-c*x+a*y)**2) +\
    (4*np.pi**3/(3*Delta**3))*((-c*eps*S+d*eps*C)*(d*x-b*y)**3 + \
                             (d*eps*S-c*Gamma)*(-c*x+a*y)**3) 
    

print('Computing derivatives')

P_x = sp.diff(P, x)
P_xx = sp.diff(P_x, x)
P_xxx = sp.diff(P_xx, x)

P_xy = sp.diff(P_x, y)
P_xyy = sp.diff(P_xy, y)

P_yy = sp.diff(sp.diff(P, y), y)

Q_x = sp.diff(Q, x)
Q_xx = sp.diff(Q_x, x)
Q_xy = sp.diff(Q_x, y)
Q_xxy = sp.diff(Q_xx, y)

Q_yy = sp.diff(sp.diff(Q, y), y)
Q_yyy = sp.diff(Q_yy, y)

lyap = (1/(8*omega))*(P_xxx + P_xyy + Q_xxy + Q_yyy) +\
       (1/(8*omega**2))*(P_xy*(P_xx+P_yy)-Q_xy*(Q_xx+Q_yy)-P_xx*Q_xx+P_yy*Q_yy)

lyap = sp.lambdify(args=[x,y,t], expr=lyap)

first_lyap = []

i = 0

for arc in t_tr0_centre:
    print(100*i/len(t_tr0_centre))
    i += 1
    
    first_lyap.append(lyap(0,0,arc))

plt.plot(t_tr0_centre, first_lyap)
plt.axhline(y=0, color='k')

plt.xlabel('$\mathbf{x}_0$')
plt.title('First Lyapunov Coefficient')  

#plt.savefig('First Lyapunov Coefficient.png', dpi = 2000)