import numpy as np
import scipy as sci

def somb(x, e):
    """SOMB ...

    SOMB 2*j1(pi*x)/(pi*x) function.
    SOMB(X) returns a matrix whose elements are the somb of the elements 
    of X, i.e.
    y = 2*j1(pi*x)/(pi*x)    if x ~= 0
    = 1                    if x == 0
    where x is an element of the input matrix and y is the resultant
    output element.  
    
    Author(s): J. Loomis, 6-29-1999
    """

    z = np.ones(x.shape, dtype=float)
    # is this finding where the values are > 0?
    i = np.where(x > 0)

    z[i] = 2.0 * (sci.special.jn(1, np.pi*x[i]) - e*scipy.special.jn(1,e*np.pi*x[i])) /  (np.pi*x[i])/(1-e**2)

    return z
