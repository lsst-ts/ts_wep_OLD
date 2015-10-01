from ZernikeAnnularFit import ZernikeAnnularFit
from ZernikeFit import ZernikeFit

import numpy as np
def ZernikeMaskedFit(S, x, y,  numTerms,mask, e ):

    j,i = np.nonzero(mask[:])
    S = S[i,j]
    x = x[i,j]
    y = y[i,j]
    if (e>0):
        Z = ZernikeAnnularFit(S, x, y, numTerms, e )
    else:
        Z = ZernikeFit(S, x, y, numTerms)
    return Z
