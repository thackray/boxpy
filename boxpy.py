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


def make_fromfile_to_fill(filename, compartments):
    """use compartment names to create a file to fill with flows"""
    return

def load_fromfile(filename):
    """Read filename to make from_i_to_j"""
    return

def find_closest(x,ylist):
    """Find the closest place in ylist to x."""
    i = np.abs(ylist-x).argmin()
    return i


def get_s_from_timeseries(taxis,timeseries,interp=True):
    """Get the model forcing function s from a given timeseries.

    Given a time axis and a values axis for those times, will
    give a functional form for use in odeint.

    Arguments:
    taxis - time axis (1D array)
    timeseries - time series values (same shape as t)
    
    Keyword Arguments:
    interp - Boolean. if True, will interpolate values. if False, will use closest value from timeseries 

    Returns:
    s - forcing function (with time as argument)
    """
    
    assert len(taxis) == len(timeseries), "Time series length must match time axis"
    assert len(set(taxis)) == len(taxis), "Time axis must not have duplicate values"

    if interp:
        if len(timeseries.shape) > 1:
            d2len = timeseries.shape[1]
        else:
            d2len = 1
        def s(tnew):
            out = np.zeros(d2len)
            for i in range(d2len):
                out[i] = np.interp(tnew,taxis,timeseries[:,i])
            return out
    else:
        def s(tnew):
            return timeseries[find_closest(tnew,taxis)]
    
    return s




if __name__=="__main__":
    
    K = np.array([[0,1.],
                  [0,-1]])
    def s(t):
        return [0,1]
    c0 = [1,0]
    t = np.linspace(0,10,10)
    solution = solve_odeint(get_f_ddt(K,s),c0,t)
    print(solution)






