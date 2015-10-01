from ZernikeEval import ZernikeEval
import numpy as np

def ZernikeFit(S, x, y, numTerms):
    #print x.shape
    #if x,y are 2D, m1,m2 are lists, still (m1!=m2) below works
    m1 = x.shape
    m2 = y.shape
    if( (m1 != m2)):
        print( 'x & y are not the same size' )

    S = S[:].copy()
    x = x[:].copy()
    y = y[:].copy()

    i = np.isfinite(S + x + y)
    S = S[i]
    x = x[i]
    y = y[i]

    H = np.zeros((len(S), int(numTerms)))

    for i in range(int(numTerms)):
        Z = np.zeros(int(numTerms))
        Z[i]=1;
        H[:,i] = ZernikeEval(Z,x,y)

    Z = np.dot(np.linalg.pinv( H ), S)

    return Z
