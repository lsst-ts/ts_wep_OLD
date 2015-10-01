import numpy
def ZernikeAnnularGrad( Z, x, y ,e, type):
    '''Gradient of the Annular Zernicke'''
    m1, n1 = x.shape
    m2, n2 = y.shape

    if( m1 != m2 or n1 != n2 ):
        print( 'x & y are not the same size' )
        exit()

    if( len(Z) > 22 ):
        print('ZernikeAnnularEval() is not implemented with >22 terms')
        return
    elif len(Z)<22:
        Z[21]=0

    x2 = x* x
    y2 = y* y
    x4 = x2*x2
    y4 = y2*y2
    xy = x* y
    r2 = x2 + y2
    r4 = r2*r2
    e2=e*e
    e4=e2*e2
    e6=e4*e2
    e8=e6*e2
    e10=e8*e2
    e12=e10*e2

    if (type== 'dx'):
        d = Z[0]  * 0 * x # to make d an array with the same size as x 
        den=numpy.sqrt(1+e2)
        d =  d + Z[1]  * 2 * 1/den
        d =  d + Z[2]  * 2 * 0
        den=1-e**2
        d =  d + Z[3]  * numpy.sqrt( 3 )  * 4 * x/den
        den=numpy.sqrt(1+e2+e4)
        d =  d + Z[4]  * numpy.sqrt( 6 )  * 2 * y/den
        d =  d + Z[5]  * numpy.sqrt( 6 )  * 2 * x/den
        den=numpy.sqrt( (1-e2)**2*(1+e2)*(1+4*e2+e4) )
        d =  d + Z[6]  * numpy.sqrt( 8 )  * 6 * xy*(1+e2)/den
        d =  d + Z[7]  * numpy.sqrt( 8 )  * ((9 * x2 + 3 * y2 - 2)*(1+e2)-2*e4)/den
        den=numpy.sqrt(1+e2+e4+e6)
        d =  d + Z[8]  * numpy.sqrt( 8 )  * 6 * xy/den
        d =  d + Z[9] * numpy.sqrt( 8 )  * (3 * x2 - 3 * y2)/den
        den=(1-e2)**2
        d =  d + Z[10] * numpy.sqrt( 5 )  * 12 * x* (2 * r2 - 1-e2)/den
        den=(1-e2)**3*(1+e2+e4)
        num=numpy.sqrt((1-e2)**4*(1+e2+e4)/(1+4*e2+10*e4+4*e6+e8))
        d =  d + Z[11] * numpy.sqrt( 10 ) * (x* (16 * x2 - 6)*(1+e2+e4)-6*x*e6)*num/den
        d =  d + Z[12] * numpy.sqrt( 10 ) * (y* (24 * x2 + 8 * y2 - 6)*(1+e2+e4)-6*y*e6)*num/den
        den=numpy.sqrt(1+e2+e4+e6+e8)
        d =  d + Z[13] * numpy.sqrt( 10 ) * 4 * x* (x2 - 3 * y2)/den
        d =  d + Z[14] * numpy.sqrt( 10 ) * 4 * y* (3 * x2 - y2)/den
        den=(1-e2)**3*(1+4*e2+e4)
        num=numpy.sqrt((1-e2)**2*(1+4*e2+e4)/(1+9*e2+9*e4+e6) )
        d =  d + Z[15] * numpy.sqrt( 12 ) * (3*e8 - 36*e6*x2 - 12*e6*y2 + 12*e6 + 50*e4*x4 + 60*e4*x2*y2 - 144*e4*x2 \
                                       + 10*e4*y4 - 48*e4*y2 + 30*e4 + 200*e2*x4 + 240*e2*x2*y2 - 144*e2*x2 + 40*e2*y4 - 48*e2*y2\
                                       + 12*e2 + 50*x4 + 60*x2*y2 - 36*x2 + 10*y4 - 12*y2 + 3 ) *num/den
        d =  d + Z[16] * numpy.sqrt( 12 ) * (8*xy* (5*r2*(1+4*e2+e4)-(3+12*e2+12*e4+3*e6) ))*num/den
        den=(1-e2)**4*(1+e2)*(1+e4)
        num=numpy.sqrt((1-e2)**6*(1+e2)*(1+e4)/(1+4*e2+10*e4+20*e6+10*e8+4*e10+e12))
        d =  d + Z[17] * numpy.sqrt( 12 ) * (25*(e6 + e4 + e2 +1)*x4 + (- 12*e8 - 30*e6*y2 - 12*e6 - 30*e4*y2 - 12*e4\
                                                                    - 30*e2*y2 - 12*e2 - 30*y2 - 12)*x2 + 12*e8*y2 - 15*e6*y4 + 12*e6*y2 - 15*e4*y4 + 12*e4*y2 \
                                       - 15*e2*y4 + 12*e2*y2 - 15*y4 + 12*y2 )*num/den
        d =  d + Z[18] * numpy.sqrt( 12 ) * (4.0*xy*(15*(e6+e4+e2+1)*x2 - 6*e8 + 5*e6*y2 - 6*e6 + 5*e4*y2 - 6*e4 + 5*e2*y2\
                                                - 6*e2 + 5*y2 - 6 ))*num/den
        den=numpy.sqrt(1+e2+e4+e6+e8+e10)
        d =  d + Z[19] * numpy.sqrt( 12 ) * 5*(x2*(x2-6*y2)+y4)/den
        d =  d + Z[20] * numpy.sqrt( 12 ) * 20*xy*(x2-y2)/den
        den=(1-e2)**3
        d =  d + Z[21] * numpy.sqrt(  7 )  * 24*x*(e4 - e2*(5*y2-3) + 5*x4 - 5*y2 + 5*y4 - x2*(5*e2 - 10*y2 + 5) + 1)/den
    elif (type == 'dy'):
        d =        Z[0]  * 0 * x
        den=numpy.sqrt(1+e2)
        d =  d + Z[1]  * 2 * 0
        d =  d + Z[2]  * 2 * 1/den
        den=1-e**2
        d =  d + Z[3]  * numpy.sqrt( 3 )  * 4 * y/den
        den=numpy.sqrt(1+e2+e4)    
        d =  d + Z[4]  * numpy.sqrt( 6 )  * 2 * x/den
        d =  d + Z[5]  * numpy.sqrt( 6 )  * (-2) * y/den
        den=numpy.sqrt( (1-e2)**2*(1+e2)*(1+4*e2+e4) )
        d =  d + Z[6]  * numpy.sqrt( 8 )  * ((1+e2)*(3 * x2 + 9 * y2 - 2) -2*e4)/den
        d =  d + Z[7]  * numpy.sqrt( 8 )  * 6 * xy*(1+e2)/den
        den=numpy.sqrt(1+e2+e4+e6)
        d =  d + Z[8]  * numpy.sqrt( 8 )  * (3 * x2 - 3 * y2)/den
        d =  d + Z[9] * numpy.sqrt( 8 )  * (-6) * xy/den
        den=(1-e2)**2
        d =  d + Z[10] * numpy.sqrt( 5 )  * 12 * y* (2 * r2 - 1-e2)/den
        den=(1-e2)**3*(1+e2+e4)
        num=numpy.sqrt((1-e2)**4*(1+e2+e4)/(1+4*e2+10*e4+4*e6+e8))
        d =  d + Z[11] * numpy.sqrt( 10 ) * (y* (6 - 16 * y2)*(1+e2+e4)+6*y*e6)*num/den
        d =  d + Z[12] * numpy.sqrt( 10 ) * (x* (8 * x2 + 24 * y2 - 6)*(1+e2+e4)-6*x*e6)*num/den
        den=numpy.sqrt(1+e2+e4+e6+e8)
        d =  d + Z[13] * numpy.sqrt( 10 ) * 4 * y* (y2 - 3 * x2)/den
        d =  d + Z[14] * numpy.sqrt( 10 ) * 4 * x* (x2 - 3 * y2)/den
        den=(1-e2)**3*(1+4*e2+e4)
        num=numpy.sqrt((1-e2)**2*(1+4*e2+e4)/(1+9*e2+9*e4+e6) )
        d =  d + Z[15] * numpy.sqrt( 12 ) * (-x*(24*y + 4*e2*(24*y - 40*y*r2) + 2*e4*(48*y - 20*y*r2) + 24*e6*y - 40*y*r2))*num/den
        d =  d + Z[16] * numpy.sqrt( 12 ) * (3*e8 - 12*e6*x2 - 36*e6*y2 + 12*e6 + 10*e4*x4 + 60*e4*x2*y2 - 48*e4*x2 \
                                       + 50*e4*y4 - 144*e4*y2 + 30*e4 + 40*e2*x4 + 240*e2*x2*y2 - 48*e2*x2 + 200*e2*y4 - 144*e2*y2 \
                                       + 12*e2 + 10*x4 + 60*x2*y2 - 12*x2 + 50*y4 - 36*y2 + 3 )*num/den
        den=(1-e2)**4*(1+e2)*(1+e4)
        num=numpy.sqrt((1-e2)**6*(1+e2)*(1+e4)/(1+4*e2+10*e4+20*e6+10*e8+4*e10+e12))    
        d =  d + Z[17] * numpy.sqrt( 12 ) * (4.0*xy*((- 5)*(e6 +e4+ e2 +1)*x2 + 6*e8 - 15*e6*y2 + 6*e6 - 15*e4*y2\
                                                + 6*e4 - 15*e2*y2 + 6*e2 - 15*y2 + 6))*num/den
        d =  d + Z[18] * numpy.sqrt( 12 ) * (- 12*e8*x2 + 12*e8*y2 + 15*e6*x4 + 30*e6*x2*y2 - 12*e6*x2 - 25*e6*y4\
                                         + 12*e6*y2 + 15*e4*x4 + 30*e4*x2*y2 - 12*e4*x2 - 25*e4*y4 + 12*e4*y2 + 15*e2*x4 + 30*e2*x2*y2\
                                         - 12*e2*x2 - 25*e2*y4 + 12*e2*y2 + 15*x4 + 30*x2*y2 - 12*x2 - 25*y4 + 12*y2 )*num/den
        den=numpy.sqrt(1+e2+e4+e6+e8+e10)
        d =  d + Z[19] * numpy.sqrt( 12 ) * 20*xy*(y2-x2)/den
        d =  d + Z[20] * numpy.sqrt( 12 ) * 5*(x2*(x2-6*y2)+y4)/den
        den=(1-e2)**3
        d =  d + Z[21] * numpy.sqrt(  7 )  * 24*y*(e4 - e2*(5*x2 - 3) - 5*x2 + 5*x4 + 5*y4 - y2*(5*e2 - 10*x2 + 5) + 1)/den
    elif (type == 'dx2'):
        d =        Z[0]  * 0 * x # to make d an array with the same size as x 
        d =  d + Z[1]  * 0
        d =  d + Z[2]  * 0
        den=1-e**2    
        d =  d + Z[3]  * numpy.sqrt( 3 )  * 4 /den
        d =  d + Z[4]  * 0
        den=numpy.sqrt(1+e2+e4)
        d =  d + Z[5]  * numpy.sqrt( 6 )  * 2 /den
        den=numpy.sqrt( (1-e2)**2*(1+e2)*(1+4*e2+e4) )
        d =  d + Z[6]  * numpy.sqrt( 8 )  * 6 * y*(1+e2)/den
        d =  d + Z[7]  * numpy.sqrt( 8 )  * 18 * x * (1+e2)/den
        den=numpy.sqrt(1+e2+e4+e6)
        d =  d + Z[8]  * numpy.sqrt( 8 )  * 6 * y/den
        d =  d + Z[9] * numpy.sqrt( 8 )  * 6 * x /den
        den=(1-e2)**2
        d =  d + Z[10] * numpy.sqrt( 5 )  * 12 * (6 * x2 + 2* y2 - e2 -1)/den
        den=(1-e2)**3*(1+e2+e4)
        num=numpy.sqrt((1-e2)**4*(1+e2+e4)/(1+4*e2+10*e4+4*e6+e8))
        d =  d + Z[11] * numpy.sqrt( 10 ) * ((48 * x2 - 6)*(1+e2+e4)-6*e6)*num/den
        d =  d + Z[12] * numpy.sqrt( 10 ) * 48 * xy *(1+e2+e4)*num/den
        den=numpy.sqrt(1+e2+e4+e6+e8)
        d =  d + Z[13] * numpy.sqrt( 10 ) * 12 *  (x2 - y2)/den
        d =  d + Z[14] * numpy.sqrt( 10 ) * 24 * xy/den
        den=(1-e2)**3*(1+4*e2+e4)
        num=numpy.sqrt((1-e2)**2*(1+4*e2+e4)/(1+9*e2+9*e4+e6) )
        d =  d + Z[15] * numpy.sqrt( 12 ) * (-8*x*(9*e6 - 25*e4*x2 - 15*e4*y2 + 36*e4 - 100*e2*x2 \
                                              - 60*e2*y2 + 36*e2 - 25*x2 - 15*y2 + 9))*num/den
        d =  d + Z[16] * numpy.sqrt( 12 ) * (-8*y*(3*e6 - 15*e4*x2 - 5*e4*y2 + 12*e4 - 60*e2*x2\
                                              - 20*e2*y2 + 12*e2 - 15*x2 - 5*y2 + 3) )*num/den
        den=(1-e2)**4*(1+e2)*(1+e4)
        num=numpy.sqrt((1-e2)**6*(1+e2)*(1+e4)/(1+4*e2+10*e4+20*e6+10*e8+4*e10+e12))
        d =  d + Z[17] * numpy.sqrt( 12 ) * ( -4*x*(6*e8 - 25*e6*x2 + 15*e6*y2 + 6*e6 - 25*e4*x2 \
                                               + 15*e4*y2 + 6*e4 - 25*e2*x2 + 15*e2*y2 + 6*e2 - 25*x2 + 15*y2 + 6))*num/den
        d =  d + Z[18] * numpy.sqrt( 12 ) * (-4*y*(6*e8 - 45*e6*x2 - 5*e6*y2 + 6*e6 - 45*e4*x2 \
                                              - 5*e4*y2 + 6*e4 - 45*e2*x2 - 5*e2*y2 + 6*e2 - 45*x2 - 5*y2 + 6))*num/den
        den=numpy.sqrt(1+e2+e4+e6+e8+e10)
        d =  d + Z[19] * numpy.sqrt( 12 ) * 20 * x* (x2-3*y2)/den
        d =  d + Z[20] * numpy.sqrt( 12 ) * 20* y*(3*x2-y2)/den
        den=(1-e2)**3
        d =  d + Z[21] * numpy.sqrt(  7 )  * (480*x2*r2 + 120*r4 + 24*e4 - 360*x2 - 120*y2 \
                                        - 3*e2*(120*x2 + 40*y2 - 24) + 24)/den
    
    elif (type == 'dy2'):
        d =        Z[0]  * 0 * x # to make d an array with the same size as x 
        d =  d + Z[1]  * 0
        d =  d + Z[2]  * 0
        den=1-e**2    
        d =  d + Z[3]  * numpy.sqrt( 3 )  * 4 /den
        d =  d + Z[4]  * 0
        den=numpy.sqrt(1+e2+e4)    
        d =  d + Z[5]  * numpy.sqrt( 6 )  * (-2) /den
        den=numpy.sqrt( (1-e2)**2*(1+e2)*(1+4*e2+e4) )
        d =  d + Z[6]  * numpy.sqrt( 8 )  * (1+e2)* 18 * y /den
        d =  d + Z[7]  * numpy.sqrt( 8 )  * 6 * x*(1+e2)/den
        den=numpy.sqrt(1+e2+e4+e6)
        d =  d + Z[8]  * numpy.sqrt( 8 )  * (-6) * y/den
        d =  d + Z[9] * numpy.sqrt( 8 )  * (-6) * x/den
        den=(1-e2)**2
        d =  d + Z[10] * numpy.sqrt( 5 )  * 12 * (2 * x2 + 6* y2 - e2 -1)/den
        den=(1-e2)**3*(1+e2+e4)
        num=numpy.sqrt((1-e2)**4*(1+e2+e4)/(1+4*e2+10*e4+4*e6+e8))
        d =  d + Z[11] * numpy.sqrt( 10 ) * ((6 - 48 * y2)*(1+e2+e4)+6*e6)*num/den
        d =  d + Z[12] * numpy.sqrt( 10 ) * 48 * xy*(1+e2+e4)*num/den
        den=numpy.sqrt(1+e2+e4+e6+e8)
        d =  d + Z[13] * numpy.sqrt( 10 ) * 12 * (y2 - x2)/den
        d =  d + Z[14] * numpy.sqrt( 10 ) * (-24) * xy/den
        den=(1-e2)**3*(1+4*e2+e4)
        num=numpy.sqrt((1-e2)**2*(1+4*e2+e4)/(1+9*e2+9*e4+e6) )
        d =  d + Z[15] * numpy.sqrt( 12 ) * (-8*x*(3*e6 - 5*e4*x2 - 15*e4*y2 + 12*e4 - 20*e2*x2\
                                              - 60*e2*y2 + 12*e2 - 5*x2 - 15*y2 + 3) ) *num/den
        d =  d + Z[16] * numpy.sqrt( 12 ) * (-8*y*(9*e6 - 15*e4*x2 - 25*e4*y2 + 36*e4 - 60*e2*x2\
                                              - 100*e2*y2 + 36*e2 - 15*x2 - 25*y2 + 9))*num/den
        den=(1-e2)**4*(1+e2)*(1+e4)
        num=numpy.sqrt((1-e2)**6*(1+e2)*(1+e4)/(1+4*e2+10*e4+20*e6+10*e8+4*e10+e12))    
        d =  d + Z[17] * numpy.sqrt( 12 ) * (4*x*(6*e8 - 5*e6*x2 - 45*e6*y2 + 6*e6 - 5*e4*x2\
                                             - 45*e4*y2 + 6*e4 - 5*e2*x2 - 45*e2*y2 + 6*e2 - 5*x2 - 45*y2 + 6))*num/den
        d =  d + Z[18] * numpy.sqrt( 12 ) * (4*y*(6*e8 + 15*e6*x2 - 25*e6*y2 + 6*e6 + 15*e4*x2 \
                                             - 25*e4*y2 + 6*e4 + 15*e2*x2 - 25*e2*y2 + 6*e2 + 15*x2 - 25*y2 + 6) )*num/den
        den=numpy.sqrt(1+e2+e4+e6+e8+e10)
        d =  d + Z[19] * numpy.sqrt( 12 ) * 20*x*(3*y2-x2)/den
        d =  d + Z[20] * numpy.sqrt( 12 ) * 20*y*(y2 - 3*x2)/den
        den=(1-e2)**3
        d =  d + Z[21] * numpy.sqrt(  7 )  * (480*y2*r2 + 120*r4 + 24*e4 - 120*x2 - 360*y2 \
                                        - 3*e2*(40*x2 + 120*y2 - 24) + 24)/den
    
    elif (type == 'dxy'):
        d =        Z[0]  * 0 * x # to make d an array with the same size as x 
        d =  d + Z[1]  * 0
        d =  d + Z[2]  * 0
        d =  d + Z[3]  * 0
        den=numpy.sqrt(1+e2+e4)
        d =  d + Z[4]  * numpy.sqrt( 6 )  * 2 /den
        d =  d + Z[5]  * 0
        den=numpy.sqrt( (1-e2)**2*(1+e2)*(1+4*e2+e4) )
        d =  d + Z[6]  * numpy.sqrt( 8 )  * (1+e2)*(6 * x)/den
        d =  d + Z[7]  * numpy.sqrt( 8 )  * 6 * y*(1+e2)/den
        den=numpy.sqrt(1+e2+e4+e6)
        d =  d + Z[8]  * numpy.sqrt( 8 )  * 6 * x/den
        d =  d + Z[9] * numpy.sqrt( 8 )  * (-6) * y/den
        den=(1-e2)**2
        d =  d + Z[10] * numpy.sqrt( 5 )  * 48 *xy/den    
        den=(1-e2)**3*(1+e2+e4)
        num=numpy.sqrt((1-e2)**4*(1+e2+e4)/(1+4*e2+10*e4+4*e6+e8))
        d =  d + Z[11] * numpy.sqrt( 10 ) * 0
        d =  d + Z[12] * numpy.sqrt( 10 ) * ((24 * x2 + 24 * y2 - 6)*(1+e2+e4)-6*e6)*num/den
        den=numpy.sqrt(1+e2+e4+e6+e8)
        d =  d + Z[13] * numpy.sqrt( 10 ) * (-24)*xy/den
        d =  d + Z[14] * numpy.sqrt( 10 ) * 12 *  (x2 - y2)/den
        den=(1-e2)**3*(1+4*e2+e4)
        num=numpy.sqrt((1-e2)**2*(1+4*e2+e4)/(1+9*e2+9*e4+e6) )
        d =  d + Z[15] * numpy.sqrt( 12 ) * (-8*y*(3*e6 - 15*e4*x2 - 5*e4*y2 + 12*e4 - 60*e2*x2\
                                              - 20*e2*y2 + 12*e2 - 15*x2 - 5*y2 + 3))*num/den
        d =  d + Z[16] * numpy.sqrt( 12 ) * (-8*x*(3*e6 - 5*e4*x2 - 15*e4*y2 + 12*e4 - 20*e2*x2\
                                              - 60*e2*y2 + 12*e2 - 5*x2 - 15*y2 + 3) )*num/den
        den=(1-e2)**4*(1+e2)*(1+e4)
        num=numpy.sqrt((1-e2)**6*(1+e2)*(1+e4)/(1+4*e2+10*e4+20*e6+10*e8+4*e10+e12))    
        d =  d + Z[17] * numpy.sqrt( 12 ) * (12*y*(2*e8 - 5*e6*r2 + 2*e6 - 5*e4*r2 + 2*e4\
                                              - 5*e2*r2 + 2*e2 - 5*r2 + 2))*num/den
        d =  d + Z[18] * numpy.sqrt( 12 ) * (-12*x*(2*e8 - 5*e6*r2 + 2*e6 - 5*e4*r2 + 2*e4 \
                                               - 5*e2*r2 + 2*e2 - 5*r2 + 2) )*num/den
        den=numpy.sqrt(1+e2+e4+e6+e8+e10)
        d =  d + Z[19] * numpy.sqrt( 12 ) * 20*y*(y2-3*x2)/den
        d =  d + Z[20] * numpy.sqrt( 12 ) * 20*x*(x2 - 3*y2)/den
        den=(1-e2)**3
        d =  d + Z[21] * numpy.sqrt(  7 )  * 240 * xy*(2*r2-1-e2)/den
    
    return d



 
