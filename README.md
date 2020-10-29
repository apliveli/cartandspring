Cart and Spring system simulations
---------------------------------------------------

About
---------

Simulation of the motion of multiple coupled-oscillator cart and spring systems using Python. This simulation produces simulated plots of the positions of coupled carts as a function of time. Carts of user defined mass, initial position, initial velocity, and friction coefficient as well as springs of user defined constants and length can be simulated.


Package requirements
--------------------------------

1. Numpy
2. matplotlib.pyplot
3. scipy.integrate

User defined variables
--------------------------------

Initial:

    Spring:
        1. First spring's constant
        2. First spring's length
    Cart:
        3. First cart's mass
        4. First cart's initial position
        5. First cart's initial velocity
        6. First cart's friction coefficient
        
After initial system is created, more springs and carts can be added using the variables listed above.

Limitations:
----------------

This simulation can handle multiple different system setups. These include:
1. One cart and one spring
2. One cart and two springs
3. Two carts and two springs
4. Two carts and three springs

Cart 2 cannot have an initial position less than cart 1. When calculating the equations of motion, cart 1 is assumed to be the first cart in the series. 
