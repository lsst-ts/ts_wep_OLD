import numpy as np

def getdIandI(m):
    
    m1, n1 = m.I1.shape
    m2, n2 = m.I2.shape

    if( m1 != n1 ):
        print( 'getdIandI: I1 is not square' )
        exit()

    if( (m1 != m2) or (n1 != n2) ):
        print( 'getdIandI: I1 and I2 are not the same size' )
        exit()

    I1 = m.I1
    I2 = np.rot90(m.I2, 2 )

    m.I = (I1+I2)/2
    m.dI = I2 - I1
    return m

