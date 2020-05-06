# Bifurcation Diagrams of BMK Equations

## Scope:
Creating the bifurcation diagram for the family of ODEs described in part 2 of Simplest Bifurcation Diagrams of Monotone Vector Fields on a Torus (C Baesens and R S MacKay, 2018). From here on, they will be called BMK equations.
This project plots the SNE curves (equivalently the boundaries of the resonance region), trace-zero loops (as continuation from bogdanov-takens points) as well as curves of RHC+/- that intersect forming the Necklace point. Ultimately, the aim is to tune the equation to form modified BMK equations and move the Necklace point inside the trace-zero loop.

## Current status
SNE and trace-zero are plotted using PyDSTool's PyCont class. 
The Ode class can find the saddle point and approximate the unstable manifold. This is computed up until it starts moving away from the (1,0)-lift of the saddle point and the distance to the lift point is computed (this distance is called the Pontryagin Energy of the system). Then apply Newton-Raphson method to find the zero (i.e. the heteroclinic connection in the plane). The argmin is a point on the RHC+ curve. 
This value is then used as the starting seed to compute the minima for an increment in ox.
The whole process is repeated from the negative branch of the unstable manifold to the (1,0)-lift of the saddle to generate RHC-.

Applying the appropriate symmetry x -> -x and y -> -y generates the RHC curves for the lower necklace point.

The N point can be computed directly by using a multivariate Newton Raphson method with a learning rate eta = 0.1. The NR method is applied to the function (ox, oy) --> (PE+, PE-). Where PE+(-) is the Pontryagin Energy for the positive (negative) branch of the unstable manifold to the ((1,0)-lift of the) saddle.


## TODO

