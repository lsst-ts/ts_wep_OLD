import numpy as np
from rotateMaskParam import rotateMaskParam
def makeCornerMasklist(obsR, ca, ra, cb, rb, fldx,fldy):

    cax, cay, cbx, cby = rotateMaskParam(ca,cb,fldx,fldy)
    masklist=np.array([[0, 0, 1, 1],[0, 0, obsR, 0],[cax, cay, ra, 1],[cbx, cby, rb, 0]])
    
    return masklist
