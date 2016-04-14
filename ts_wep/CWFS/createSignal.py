import numpy as np
from padArray import padArray

def createSignal( m , cliplevel):

    # missing first line in the image

    m1, n1 = m.I1.shape
    m2, n2 = m.I2.shape

    if( m1 != n1 ):
        raise Exception( 'EFSignal: I1 is not square' )

    if( (m1 != m2) or (n1 != n2) ):
        raise Exception( 'EFSignal: I1 and I2 are not the same size' )

    I1 = m.I1
    I2 = np.rot90(m.I2.copy(), k=2 ) #do not change m.I2 in PoissionSolver.m (

    num = I1 - I2 # -(I2-I1), the - is from S itself, see Eq.(4) of our SPIE
    den = I1 + I2

    #to apply signal_sum_clip_level
    pixelList=den*m.cMask
    pixelList[pixelList==0]=np.nan
    m1,n1= m.cMask.shape
    pixelList = np.reshape(pixelList,m1*n1)
    pixelList = pixelList[~np.isnan(pixelList)]
    low= pixelList.min()
    high= pixelList.max()
    median = (high-low)/2.+low
    den[den<median*cliplevel] = 1.5*median

    i =  den[:] == 0
    den[i] = np.inf # Forces zero in the result below.
    m.S = num / den

    c0=m.focalLength*(m.focalLength-m.offset)/m.offset
    m.S = m.S/c0

    m.S = padArray( m.S, m.padDim )*m.cMaskPad

    return m
