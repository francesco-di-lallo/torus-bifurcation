# Bifurcation Diagrams of BMK Equations

## Scope:
Creating the bifurcation diagram for the family of ODEs described in part 2 of Simplest Bifurcation Diagrams of Monotone Vector Fields on a Torus (C Baesens and R S MacKay, 2018). From here on, they will be called BMK equations.
This project plots the SNE curves (equivalently the boundaries of the resonance region), trace zero loops as continuation from bogdanov-takens points as well as curves of RHC+/- that intersect forming the Necklace point. Ultimately, the aim is to tune the equation to form modified BMK equations and move the Necklace point inside the tr0 loop.

## Current status
SNE and Tr0 are plotted using PyDSTool's PyCont class. 
The Ode class can find the saddle point and approximate the unstable manifold. This is computed up until it starts moving away from the (1,0)-lift of the saddle point and the distance to the lift point is computed (this distance is called the Pontryagin Energy of the system). Then apply a gradient descent algorithm to find oy that minimises PE. The argmin is a point on the RHC curve. 
This value is then used as the starting seed to compute the minima for an increment in ox.

## Progression
Repeating this whole process but initialising from the negative branch unstable manifold of the (1,0)-lift of the saddle point will generate the RHC- curve.

## TODO
- Correct gradient descent algorithm since PE curve is not differentiable at minimum --> Try transforming the distance into a smoother function. Will look into squaring the list or exponentiating it.
- Propagate the unstable manifold of the (1,0)-lift of the saddle backwards that coincides with the forward branch of the stable manifold the saddle point at a RHC point. This generates the RHC- curve.
