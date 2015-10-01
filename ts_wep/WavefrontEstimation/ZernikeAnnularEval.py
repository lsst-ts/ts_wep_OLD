import numpy as np 

def ZernikeAnnularEval(Z, x, y, e):
    '''Evaluate the Annular Zernickes'''

    m1 = x.shape
    m2 = y.shape

    if (m1 != m2):
        print( 'x & y are not the same size' )
        exit()

    if( len(Z) > 22 ):
        print('ZernikeAnnularEval() is not implemented with >22 terms')
        return;
    elif len(Z)<22:
        Z[21]=0

    r2 = x*x + y*y
    r = np.sqrt( r2 )
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

    e2=e*e
    e4=e2*e2
    e6=e4*e2
    e8=e6*e2
    e10=e8*e2
    e12=e10*e2

    S = Z[0] * (1 + 0*x) # 0*x to set NaNs properly

    den=np.sqrt(1+e2)
    S = S + Z[1]  * 2 * r * c/den
    S = S + Z[2]  * 2 * r * s/den

    den=1-e**2
    S = S + Z[3]  * np.sqrt( 3 )  * (2*r2 - 1-e2)/den

    den=np.sqrt(1+e2+e4)
    S = S + Z[4]  * np.sqrt( 6 )  * r2 * s2/den
    S = S + Z[5]  * np.sqrt( 6 )  * r2 * c2/den

    den=np.sqrt( (1-e2)**2*(1+e2)*(1+4*e2+e4) )
    S = S + Z[6]  * np.sqrt( 8 )  * (3*r3 - 2*r -2*e4*r + e2*r*(3*r2-2) ) * s/den
    S = S + Z[7]  * np.sqrt( 8 )  * (3*r3 - 2*r -2*e4*r + e2*r*(3*r2-2) ) * c/den

    den=np.sqrt(1+e2+e4+e6)
    S = S + Z[8]  * np.sqrt( 8 )  * r3* s3/den
    S = S + Z[9] * np.sqrt( 8 )  * r3* c3/den

    den=(1-e2)**2
    S = S + Z[10] * np.sqrt( 5 )  * (6*r4 - 6*r2 + 1+e4+e2*(4-6*r2) )/den

    den=(1-e2)**3*(1+e2+e4)
    num=np.sqrt((1-e2)**4*(1+e2+e4)/(1+4*e2+10*e4+4*e6+e8))
    S = S + Z[11] * np.sqrt( 10 ) * (4*r4 - 3*r2-3*e6*r2-e2*r2*(3-4*r2)-e4*r2*(3-4*r2)) * c2*num/den
    S = S + Z[12] * np.sqrt( 10 ) * (4*r4 - 3*r2-3*e6*r2-e2*r2*(3-4*r2)-e4*r2*(3-4*r2)) * s2*num/den

    den=np.sqrt(1+e2+e4+e6+e8)
    S = S + Z[13] * np.sqrt( 10 ) * r4 * c4/den
    S = S + Z[14] * np.sqrt( 10 ) * r4 * s4/den

    den=(1-e2)**3*(1+4*e2+e4)
    numE=np.sqrt((1-e2)**2*(1+4*e2+e4)/(1+9*e2+9*e4+e6) )
    numQ=10*r5-12*r3 +3*r+3*e8*r-12*e6*r*(r2-1)+2*e4*r*(15-24*r2+5*r4)+4*e2*r*(3-12*r2+10*r4)
    S = S + Z[15] * np.sqrt( 12 ) * numE*numQ* c/den
    S = S + Z[16] * np.sqrt( 12 ) * numE*numQ* s/den

    den=(1-e2)**4*(1+e2)*(1+e4)
    numE=np.sqrt((1-e2)**6*(1+e2)*(1+e4)/(1+4*e2+10*e4+20*e6+10*e8+4*e10+e12))
    numQ=r3*(5*r2-4-4*e8-e2*(4-5*r2)-e4*(4-5*r2)-e6*(4-5*r2))
    S = S + Z[17] * np.sqrt( 12 ) * numE*numQ* c3/den
    S = S + Z[18] * np.sqrt( 12 ) * numE*numQ* s3/den

    den=np.sqrt(1+e2+e4+e6+e8+e10)
    S = S + Z[19] * np.sqrt( 12 ) * r5 * c5/den
    S = S + Z[20] * np.sqrt( 12 ) * r5 * s5/den

    den=(1-e2)**3
    S = S + Z[21] * np.sqrt(  7 ) * (20*r6 - 30*r4 + 12*r2 - 1-e6+3*e4*(-3+4*r2)-3*e2*(3-12*r2+10*r4))/den

    return S
