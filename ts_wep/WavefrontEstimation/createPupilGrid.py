import numpy as np
from rotateMaskParam import rotateMaskParam

def  createPupilGrid(lutx,luty,onepixel,ca,cb,ra,rb,fldx=None,fldy=None):
    """Create the pupil grid"""

    if ((fldx == None) and (fldy ==None)):
        fldr=fldx
        fldx=fldr/1.4142
        fldy=fldx

    cax, cay, cbx, cby = rotateMaskParam(ca,cb,fldx,fldy)

    lutr = np.sqrt((lutx-cax)**2 + (luty-cay)**2)
    tmp=lutr.copy()
    tmp[np.isnan(tmp)]=-999
    idxout = (tmp>ra+onepixel) 
    lutx[idxout] = np.nan
    luty[idxout] = np.nan
    idxbound =  (tmp<=ra+onepixel) & (tmp>ra) & (~np.isnan(lutx)) 
    lutx[idxbound] = (lutx[idxbound]-cax)/lutr[idxbound]*ra+cax
    luty[idxbound] = (luty[idxbound]-cay)/lutr[idxbound]*ra+cay

    lutr=np.sqrt((lutx-cbx)**2 + (luty-cby)**2)
    tmp=lutr.copy()
    tmp[np.isnan(tmp)]=999
    idxout =  (tmp<rb-onepixel)
    lutx[idxout] = np.nan
    luty[idxout] = np.nan
    idxbound =   (tmp>=rb-onepixel) & (tmp<rb) & (~np.isnan(lutx))
    lutx[idxbound] = (lutx[idxbound]-cbx)/lutr[idxbound]*rb+cbx
    luty[idxbound] = (luty[idxbound]-cby)/lutr[idxbound]*rb+cby
    

    return lutx, luty



