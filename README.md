# Bifurcation Diagrams of BMK Equations

## Scope:
Creating the bifurcation diagram for the family of ODEs described in part 2 of Simplest Bifurcation Diagrams of Monotone Vector Fields on a Torus (C Baesens and R S MacKay, 2018). 
This project plots the SNE curves (equivalently the boundaries of the resonance region), trace zero loops as continuation from bogdanov-takens points as well as curves of RHC+/-.

## Current status
SNE and Tr0 are plotted using PyDSTool's PyCont class. 
The Ode class can find the saddle point and approximate the unstable manifold. This is computed up until it starts moving away from the lift of the saddle point and the distance to the lift point is computed.

## Progression
The next step is to vary the parameters of the ODE and recomputing the distance until a point with heteroclinic connection is found. Then continuing the curve to form either RHC+ or RHC-. The finding a point on the other curve and continuing that to find the N point.
Ultimately, we will modify the BMK equation to move the N point inside the Tr0 loop.

