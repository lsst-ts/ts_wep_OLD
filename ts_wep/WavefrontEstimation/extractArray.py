import numpy as np
def extractArray(inArray, dim ):

    m, n = inArray.shape
    if m != n:
        print( 'extractArray: array is not square' )

    if m < dim:
        print( 'extractArray: array is smaller than dimension' );

    #print "DIMEN", dim
    i = np.floor((m - dim)/2)
    j = i + dim;
    out = inArray[i:j, i:j]

    
    return out
