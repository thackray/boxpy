""" Set up and solve the multi-compartment box model equation
d(c)/dt = K c + s(t)
for arbitrary len(c), K, s(t)
"""
from __future__ import print_function
from scipy.integrate import odeint
import numpy as np

def solve_odeint(f_ddt, c0, t, ):
    """Solve the box-model equation using scipy's odeint.
    
    Arguments:
    f_ddt - function which computes the time derivative
    c0 - initial conditions for each compartment
    t - timeaxis of desired solution

    Returns:
    c(t)
    """

    c = odeint(f_ddt, c0, t)

    return c


def get_f_ddt(K, s):
    """Generate the f_ddt function to pass to the odeint solver.

    Assuming ddt form: dc/dt = Kc +s
    for rates of change in vector c of compartments 

    Arguments:
    K - matrix of flows; shape = (len(c),len(c))
    s - vector of forcings; shape = len(c)

    Returns:
    f_ddt - function to be passed to odeint
    """

    def f_ddt(c,t):
        dcdt = np.dot(K,c)+s
        return dcdt

    return f_ddt


if __name__=="__main__":
    
    K = np.array([[-1.1,0.],
                  [0.1,0]])
    s = [1,0]
    c0 = [1,0]
    t = np.linspace(0,10,10)
    solution = solve_odeint(get_f_ddt(K,s),c0,t)
    print(solution)






