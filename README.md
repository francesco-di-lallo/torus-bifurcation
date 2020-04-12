# Bifurcation Diagrams of BMK Equations

## Scope:
Creating the bifurcation diagram for the family of ODEs described in part 2 of Simplest Bifurcation Diagrams of Monotone Vector Fields on a Torus (C Baesens and R S MacKay, 2018). From here on, they will be called BMK equations.
This project plots the SNE curves (equivalently the boundaries of the resonance region), trace zero loops as continuation from bogdanov-takens points as well as curves of RHC+/-.

## Current status
SNE and Tr0 are plotted using PyDSTool's PyCont class. 
The Ode class can find the saddle point and approximate the unstable manifold. This is computed up until it starts moving away from the lift of the saddle point and the distance to the lift point is computed.

## Progression
The next step is to vary the parameters of the BMK class and recomputing the distance until a point with heteroclinic connection is found. This implies that under the torus mapping there is a rotational homoclinic connection. Then continuing the curve to form either RHC+ or RHC-. The finding a point on the other curve and continuing that to find the N point.
Ultimately, we will modify the BMK equation to move the N point inside the Tr0 loop.

## Technical Details
### `bmk.py` 
This is the master file that creates the bifurcation diagram. 
We initialise the ODE with given arguements. 
### Functions:
- `inner_outer_init(type PyDSTool.ContClass)`
> Computes two points: one on both boundaries of SNE curves.

- `inner_sne(type PyDSTool.ContClass)`
> Continues the inner boundary of resonance region as the continuation of a saddle-node bifurcation. Detects and labels Bogdanov-Takens points.
- `outer_sne(type PyDSTool.ContClass)`
> Continues the outer boundary of resonance region as the continuation of a saddle-node bifurcation. Detects and labels Bogdanov-Takens points.
- `trace_zero_upper(type PyDSTool.ContClass)`
> Creates upper curve of trace 0 emanating from Bogdanov-Takens points on the SNE curve. Continuation of Neutral Saddle type.
- `trace_zero_lower(type PyDSTool.ContClass)`
> Creates lower curve of trace 0 emanating from Bogdanov-Takens points on the SNE curve. Continuation of Neutral Saddle type.

Due to point inheritance, these functions have to be executed with respect to the following hierarchy:  
`inner_outer_init`  
+-- `inner_sne`  
|   +-- `trace_zero_upper`  
|   +-- `trace_zero_lower`  
+-- `inner_sne`  

Displaying the bifurcation diagram is then done by executing `PC.display(['ox','oy'])`

### `rhc_cont.py` 
In this script, we define two classes: `Ode(object)` and `RhcCont(PyDSTool.ContClass)`.  
- `Ode(object)`
Methods:  

\_\_init\_\_(ode)  

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
\_dist(x,y)  

        dist is the L2 distance between x and y

        Parameters
        ----------
        x : point dict
            
        y : point dict
            

        Returns
        -------
        norm : float
            L2 norm of x, y.

\_initialise_unstable_manifold(perturbation = 1e-5)

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
         
\_lift(fp)  

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

\_point_dict_getter(point_dict)  

        Converts point dictionary to a list

        Parameters
        ----------
        point_dict : dict

        Returns
        -------
        list
\_point_dict_setter(x,y)

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

\_RK4\_iter(point_dict, h=1e-0)

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
        dict
            Next iterate of RK4. 

evaluate(x,y)

        Evaluates the ODE as a vector field.

        Parameters
        ----------
        x : float

        y : float

        Returns
        -------
        array
## TODO
- Make a class that initiates a BMK with given parameters.
