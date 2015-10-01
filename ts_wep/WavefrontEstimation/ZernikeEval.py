import numpy as np
def ZernikeEval( Z, x, y ):
    '''Evaluate Zernicke'''

    #if x,y are 2D, m1,m2 are lists, still (m1!=m2) below works
    m1 = x.shape
    m2 = y.shape

    #print Z.shape
    #print Z
    
    if( (m1 != m2)):
        print( 'x & y are not the same size' )
        exit()

    if( len(Z) > 22 ):
        print('ZernikeAnnularEval() is not implemented with >22 terms')
        return
    elif len(Z)<22:
        Z[21]=0

    r2 = x*x + y*y
    r = np.sqrt(r2)
    r3 = r2 * r
    r4 = r2 * r2
    r5 = r3 * r2
    r6 = r3 * r3

    t = np.arctan2( y, x )
    s = np.sin( t )
    c = np.cos( t )
    s2 = np.sin( 2*t )
    c2 = np.cos( 2*t )
    s3 = np.sin( 3*t )
    c3 = np.cos( 3*t )
    s4 = np.sin( 4*t )
    c4 = np.cos( 4*t )
    s5 = np.sin( 5*t )
    c5 = np.cos( 5*t )

    S =     Z[0]  * (1 + 0*x) # 0*x to set NaNs properly
    S = S + Z[1]  * 2 * r * c
    S = S + Z[2]  * 2 * r * s
    S = S + Z[3]  * np.sqrt( 3 )  * (2*r2 - 1)
    S = S + Z[4]  * np.sqrt( 6 )  * r2 * s2
    S = S + Z[5]  * np.sqrt( 6 )  * r2 * c2
    S = S + Z[6]  * np.sqrt( 8 )  * (3*r3 - 2*r) * s
    S = S + Z[7]  * np.sqrt( 8 )  * (3*r3 - 2*r) * c
    S = S + Z[8]  * np.sqrt( 8 )  * r3* s3
    S = S + Z[9] * np.sqrt( 8 )  * r3* c3
    S = S + Z[10] * np.sqrt( 5 )  * (6*r4 - 6*r2 + 1)
    S = S + Z[11] * np.sqrt( 10 ) * (4*r4 - 3*r2) * c2
    S = S + Z[12] * np.sqrt( 10 ) * (4*r4 - 3*r2) * s2
    S = S + Z[13] * np.sqrt( 10 ) * r4 * c4
    S = S + Z[14] * np.sqrt( 10 ) * r4 * s4
    S = S + Z[15] * np.sqrt( 12 ) * (10*r5 - 12*r3 + 3*r) * c
    S = S + Z[16] * np.sqrt( 12 ) * (10*r5 - 12*r3 + 3*r) * s
    S = S + Z[17] * np.sqrt( 12 ) * (5*r5 - 4*r3) * c3
    S = S + Z[18] * np.sqrt( 12 ) * (5*r5 - 4*r3) * s3
    S = S + Z[19] * np.sqrt( 12 ) * r5 * c5
    S = S + Z[20] * np.sqrt( 12 ) * r5 * s5
    S = S + Z[21] * np.sqrt(  7 ) * (20*r6 - 30*r4 + 12*r2 - 1)


    return S
