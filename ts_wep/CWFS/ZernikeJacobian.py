import numpy as np
def ZernikeJacobian(Z, x, y , atype):

    m1, n1 = x.shape
    m2, n2 = y.shape
    if( m1 != m2 or n1 != n2 ):
        print( 'x & y are not the same size' )

    if( len(Z) > 22 ) :
        disp('ZernikeGrad() is not implemented with >22 terms')
        return
    elif len(Z)<22:
        Z[21]=0

    x2 = x* x
    y2 = y* y
    xy = x* y
    r2 = x2 + y2


    if (atype == '1st'):
        j =        Z[0]  * 0 * x # to make d an array with the same size as x
        j =  j + Z[1] * 0
        j =  j + Z[2]  * 0
        j =  j + Z[3]  * np.sqrt( 3 )  * 8
        j =  j + Z[4] * np.sqrt( 6 )  * 0
        j =  j + Z[5] * np.sqrt( 6 )  * 0
        j =  j + Z[6] * np.sqrt( 8 )  * 24 * y #W8 in Roddier's 1993 table
        j =  j + Z[7] * np.sqrt( 8 )  * 24 * x #W7 in Roddier's 1993 table
        j =  j + Z[8] * np.sqrt( 8 )  * 0
        j =  j + Z[9] * np.sqrt( 8 )  * 0
        j =  j + Z[10] * np.sqrt( 5 )  * (96 * r2 - 24)
        j =  j + Z[11] * np.sqrt( 10 ) * 48 * (x2 - y2) 
        j =  j + Z[12] * np.sqrt( 10 ) * 96 * xy
        j =  j + Z[13] * np.sqrt( 10 ) * 0
        j =  j + Z[14] * np.sqrt( 10 ) * 0
        j =  j + Z[15] * np.sqrt( 12 ) * x*(240.0*r2-96.0)
        j =  j + Z[16] * np.sqrt( 12 ) * y*(240.0*(x2+y2)-96.0)
        j =  j + Z[17] * np.sqrt( 12 ) * 80.0*x*(x2-3.0*y2)
        j =  j + Z[18] * np.sqrt( 12 ) * 80.0*y*(3*x2-y2)
        j =  j + Z[19] * np.sqrt( 12 ) * 0
        j =  j + Z[20] * np.sqrt( 12 ) * 0
        j =  j + Z[21] * np.sqrt(  7 )  * 48*(1+x2*(30*y2+15*x2-10)+y2*(15*y2-10))


    elif (atype == '2nd'):
    
        j =        Z[0]**2  * 0 * x # to make d an array with the same size as x
        j =  j + Z[1]**2  * 0
        j =  j + Z[2]**2  * 0
        j =  j + Z[3]**2 * ( 3 )  * 16
        j =  j + Z[4]**2 * ( 6 )  * (-4)
        j =  j + Z[5]**2 * ( 6 )  * (-4)
        j =  j + Z[6]**2 * ( 8 )  * (108 * y2 - 36 * x2) #W8 in Roddier's 1993 table
        j =  j + Z[7]**2 * ( 8 )  * (108 * x2 - 36 * y2) #W7 in Roddier's 1993 table
        j =  j + Z[8]**2 * ( 8 )  * (-36 * r2)
        j =  j + Z[9]**2 * ( 8 )  * (-36 * r2)
        j =  j + Z[10]**2 * ( 5 )  * 144 * (12 * r2**2 - 8 * r2 + 1)
        j =  j + Z[11]**2 * ( 10 ) * 36 * (8 * x2 - 1)* (1 - 8 * y2)
        j =  j + Z[12]**2 * ( 10 ) * 36 * (8 * r2 - 16 * (x2 - y2)**2 - 1)
        j =  j + Z[13]**2 * ( 10 ) * (-144) * r2**2
        j =  j + Z[14]**2 * ( 10 ) * (-144) * r2**2
        j =  j + Z[15]**2 * ( 12 ) * 64*(5.0*(x2+y2)-3)*(x2*(25.0*x2+20.0*y2-9)-y2*(5.0*y2-3.0))
        j =  j + Z[16]**2 * ( 12 ) * 64*(5.0*(x2+y2)-3)*(y2*(25.0*y2+20.0*x2-9)-x2*(5.0*x2-3.0))
        j =  j + Z[17]**2 * ( 12 ) * 16.0*(x2*(-36.0+360*y2+x2*(180-1275*y2-125*x2))+y2*(-36+y2*(180+225*(x2-y2))))
        j =  j + Z[18]**2 * ( 12 ) * 16.0*(y2*(-36.0+360*x2+y2*(180-1275*x2-125*y2))+x2*(-36+x2*(180+225*(y2-x2))))
        j =  j + Z[19]**2 * ( 12 ) * (-400) * r2**3
        j =  j + Z[20]**2 * ( 12 ) * (-400) * r2**3
        j =  j + Z[21]**2 * (  7 )  * 576*(1+x2*(10*y2+5*x2-5)+y2*(5*y2-5))*(1+x2*(25*x2+50*y2-15)+y2*(25*y2-15))
        

    return j
