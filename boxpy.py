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
    s - function to return vector of forcings; shape = len(c)

    Returns:
    f_ddt - function to be passed to odeint
    """

    def f_ddt(c,t):
        dcdt = np.dot(K,c)+s(t)
        return dcdt

    return f_ddt

def construct_K(from_i_to_j, names):
    """Construct the model operator K from dictionary of flows.

    Use the dictionary from_i_to_j and the list of compartment names
    to construct the K operator necessary.

    Arguments:
    from_i_to_j - dictionary with the structure { 'comp1' : {'comp1:X1,...'compN:XN}, ... 'compN':{...}} where the upper level dictionary is where the fluxes are from and the sub-dictionaries are where the flows are to
    names - list of the compartment names used in from_i_to_j

    Returns:
    K - matrix operator for f_ddt
    """
    K = np.zeros((len(names),len(names)))
    for i,namei in enumerate(names):
        for j,namej in enumerate(names):
            K[j,i] += from_i_to_j[namei][namej]
            if i != j:
                K[i,i] -= from_i_to_j[namei][namej]
    return K


if __name__=="__main__":
    
    K = np.array([[0,1.],
                  [0,-1]])
    def s(t):
        return [0,1]
    c0 = [1,0]
    t = np.linspace(0,10,10)
    solution = solve_odeint(get_f_ddt(K,s),c0,t)
    print(solution)






