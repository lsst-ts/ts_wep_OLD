import numpy as np
def rotateMaskParam(ca,cb,fldx,fldy):

    fldr=np.sqrt(fldx**2+fldy**2)
    c = fldx/fldr
    s = fldy/fldr
    
    cax = c*ca
    cay = s*ca
    cbx = c*cb
    cby = s*cb

    return cax, cay, cbx, cby
