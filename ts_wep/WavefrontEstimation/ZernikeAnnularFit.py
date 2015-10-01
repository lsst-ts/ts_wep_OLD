from ZernikeAnnularEval import ZernikeAnnularEval
import numpy as np
def ZernikeAnnularFit( S, x, y,numTerms ,e):

    m1 = x.shape
    m2 = y.shape
    if (m1 != m2):
        print( 'x & y are not the same size' )

    S = S[:].copy()
    x = x[:].copy()
    y = y[:].copy()

    i = np.isfinite(S + x + y)
    S = S[i]
    x = x[i]
    y = y[i]

    H = np.zeros((len(S), numTerms))

    for i in range(numTerms):
        Z = np.zeros((numTerms))
        Z[i]=1;
        H[:,i] = ZernikeAnnularEval(Z,x,y,e)

    Z = np.dot(np.linalg.pinv( H ), S)

    return Z
