import numpy as np
def padArray( inArray, dim ):
    m, n = inArray.shape
    if (m != n):
        raise( 'padArray: array is not square' )

    if m > dim:
        raise( 'padArray: array is larger than dimension' );

    out = np.zeros( (dim, dim) )
    i = np.floor((dim - m)/2)
    j = i + m 
    out[i:j, i:j] = inArray;

    return out
