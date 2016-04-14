import numpy as np
def ZernikeGrad( Z, x, y , atype):

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

    if (atype == 'dx'):
        d =      Z[0]  * 0 * x # to make d an array with the same size as x
        d =  d + Z[1]  * 2 * 1
        d =  d + Z[2]  * 2 * 0
        d =  d + Z[3]  * np.sqrt( 3 )  * 4 * x
        d =  d + Z[4]  * np.sqrt( 6 )  * 2 * y
        d =  d + Z[5]  * np.sqrt( 6 )  * 2 * x
        d =  d + Z[6]  * np.sqrt( 8 )  * 6 * xy
        d =  d + Z[7]  * np.sqrt( 8 )  * (9 * x2 + 3 * y2 - 2)
        d =  d + Z[8]  * np.sqrt( 8 )  * 6 * xy
        d =  d + Z[9] * np.sqrt( 8 )  * (3 * x2 - 3 * y2)
        d =  d + Z[10] * np.sqrt( 5 )  * 12 * x* (2 * (x2 + y2) - 1)
        d =  d + Z[11] * np.sqrt( 10 ) * x* (16 * x2 - 6)
        d =  d + Z[12] * np.sqrt( 10 ) * y* (24 * x2 + 8 * y2 - 6)
        d =  d + Z[13] * np.sqrt( 10 ) * 4 * x* (x2 - 3 * y2)
        d =  d + Z[14] * np.sqrt( 10 ) * 4 * y* (3 * x2 - y2)
        d =  d + Z[15] * np.sqrt( 12 ) * (x2* (50.0*x2+60.0*y2-36.0)+y2*(10.0*y2-12.0)+3)
        d =  d + Z[16] * np.sqrt( 12 ) * (xy* (40.0*r2-24.0))
        d =  d + Z[17] * np.sqrt( 12 ) * (x2* (25.0*x2-12.0-30.0*y2)+y2*(12.0-15.0*y2))
        d =  d + Z[18] * np.sqrt( 12 ) * (4.0*xy*(-6.0+15.0*x2+5.0*y2))
        d =  d + Z[19] * np.sqrt( 12 ) * 5*(x2*(x2-6*y2)+y2*y2)
        d =  d + Z[20] * np.sqrt( 12 ) * 20*xy*(x2-y2)
        d =  d + Z[21] * np.sqrt(  7 )  * 24*x*(1+x2*(10*y2-5+5*x2)+y2*(5*y2-5))
        
    elif (atype, 'dy'):
    
        d =      Z[0]  * 0 * x
        d =  d + Z[1]  * 2 * 0
        d =  d + Z[2]  * 2 * 1
        d =  d + Z[3]  * np.sqrt( 3 )  * 4 * y
        d =  d + Z[4]  * np.sqrt( 6 )  * 2 * x
        d =  d + Z[5]  * np.sqrt( 6 )  * (-2) * y
        d =  d + Z[6]  * np.sqrt( 8 )  * (3 * x2 + 9 * y2 - 2)
        d =  d + Z[7]  * np.sqrt( 8 )  * 6 * xy
        d =  d + Z[8]  * np.sqrt( 8 )  * (3 * x2 - 3 * y2)
        d =  d + Z[9] * np.sqrt( 8 )  * (-6) * xy
        d =  d + Z[10] * np.sqrt( 5 )  * 12 * y* (2 * (x2 + y2) - 1)
        d =  d + Z[11] * np.sqrt( 10 ) * y* (6 - 16 * y2)
        d =  d + Z[12] * np.sqrt( 10 ) * x* (8 * x2 + 24 * y2 - 6)
        d =  d + Z[13] * np.sqrt( 10 ) * 4 * y* (y2 - 3 * x2)
        d =  d + Z[14] * np.sqrt( 10 ) * 4 * x* (x2 - 3 * y2)
        d =  d + Z[15] * np.sqrt( 12 ) * (xy* (40.0*r2-24.0))
        d =  d + Z[16] * np.sqrt( 12 ) * (x2*(10.0*x2+60.0*y2-12.0)+y2*(50.0*y2-36.0)+3)
        d =  d + Z[17] * np.sqrt( 12 ) * (4.0*xy*(6.0-5.0*x2-15.0*y2))
        d =  d + Z[18] * np.sqrt( 12 ) * (y2*(-25.0*y2+12.0+30.0*x2)+x2*(-12.0+15.0*x2))
        d =  d + Z[19] * np.sqrt( 12 ) * 20*xy*(y2-x2)
        d =  d + Z[20] * np.sqrt( 12 ) * 5*(x2*(x2-6*y2)+y2*y2)
        d =  d + Z[21] * np.sqrt(  7 )  * 24*y*(1+y2*(10*x2-5+5*y2)+x2*(5*x2-5))
        

# 	int grad(const int n,const float x,const float y,float *dx,float *dy)
# 	     //
# 	     // return the zernike gradient,
# 	     // even terms are cosine (x), odd terms are sine (y).
# 	     //
# 	{
# 		int err = 0;
# 		*dx = *dy = 0.0;
# 		switch (n) {
# 			case (1): // piston
# 				break;
# 			case (2): // x-tilt
# 				*dx = 1.0;
# 				break;
# 			case (3): // y-tilt
# 				*dy = 1.0;
# 				break;
# 			case (4): // defocus
# 				*dx = 4.0*x;
# 				*dy = 4.0*y;
# 				break;
# 			case (5): // y-astigmatism
# 				*dx = 2.0*y;
# 				*dy = 2.0*x;
# 				break;
# 			case (6): // x-astigmatism
# 				*dx = 2.0*x;
# 				*dy = -2.0*y;
# 				break;
# 			case (7): // y-coma
# 				*dx = 6.0*x*y;
# 				*dy = 3*x*x+9*y*y-2.0;
# 				break;
# 			case (8): // x-coma
# 				*dx = 9*x*x+3*y*y-2.0;
# 				*dy = 6.0*x*y;
# 				break;
# 			case (9): // y-tri-coma
# 				*dx = 6.0*x*y;
# 				*dy = 3*x*x-3*y*y;
# 				break;
# 			case (10): // x-tri-coma
# 				*dx = 3*x*x-3*y*y;
# 				*dy = -6.0*x*y;
# 				break;
# 			case (11):{ // spherical
# 				float r2 = 2.0*(x*x+y*y)-1.0;
# 				*dx = 12.0*x*r2;
# 				*dy = 12.0*y*r2;
# 			} break;
# 			case (12): //x-5th order astigmatism
# 				*dx = x*(16*x*x-6);
# 				*dy = y*(6-16*y*y);
# 				break;
# 			case (13): //y-5th order astigmatism
# 				*dx = y*(24*x*x+8*y*y-6);
# 				*dy = x*(8*x*x+24*y*y-6);
# 				break;
# 			case (14): //x-quad  astigmatism
# 				*dx = 4*x*(x*x-3*y*y);
# 				*dy = 4*y*(y*y-3*x*x);
# 				break;
# 			case (15): //y-quad  astigmatism
# 				*dx = 4*y*(3*x*x-y*y);
# 				*dy = 4*x*(x*x-3*y*y);
# 				break;
# 			case (16):{
# 				float x2 = x*x,y2 = y*y;
# 				*dx = x2*(50.0*x2+60.0*y2-36.0)+y2*(10.0*y2-12.0)+3;
# 				*dy = x*y*(40.0*y2+40.0*x2-24.0);
# 			}break;
# 			case (17):{
# 				float x2 = x*x,y2 = y*y;
# 				*dx = x*y*(40.0*x2+20.0*y2-24.0);
# 				*dy = x2*(10.0*x2+60.0*y2-12.0)+y2*(50.0*y2-36.0)+3;
# 			}break;
# 			case (18):{
# 				float x2 = x*x,y2 = y*y;
# 				*dx = x2*(25.0*x2-12.0-30.0*y2)+y2*(12.0-15.0*y2);
# 				*dy = 4.0*x*y*(6.0-5.0*x2-15.0*y2);
# 			}break;
# 			case (19):{
# 				float x2 = x*x,y2 = y*y;
# 				*dx = 4.0*x*y*(-6.0+5.0*x2+15.0*y2);
# 				*dy = y2*(-25.0*y2+12.0+30.0*x2)+x2*(-12.0-15.0*x2);
# 			}break;
# 			case (20):{
# 				float x2 = x*x,y2 = y*y;
# 				*dx = 5*(x2*(x2-6*y2)+y2*y2);
# 				*dy = 20*x*y*(y2-x2);
# 			}break;
# 			case (21):{
# 				float x2 = x*x,y2 = y*y;
# 				*dx = 20*x*y*(x2-y2);
# 				*dy = 5*(x2*(x2-6*y2)+y2*y2);
# 			}break;
# 			case (22):{
# 				float x2 = x*x,y2 = y*y;
# 				*dx = 24*x*(1+x2*(10*y2-5+5*x2)+y2*(5*y2-5));
# 				*dy = 24*y*(1+y2*(10*x2-5+5*y2)+x2*(5*x2-5));
# 			}break;
# 			case (23):{ // y n^th order astigmatism
# 			float x2 = x*x,y2 = y*y;
# 				*dx = 2*y*(6+y2*(-20+15*y2)+x2*(-60 +75*x2+90*y2));
# 				*dy = 2*x*(6+x2*(-20+15*x2)+y2*(-60+75*y2+90*x2));
# 			}break;
# 			case (24):{ // x n^th order astigmatism
# 				float x2 = x*x,y2 = y*y;
# 				*dx = 2*x*(6+x2*(-40.0+45*x2+30*y2)-15*y2*y2);
# 				*dy = 2*y*(-6.0+y2*(40-30*x2-45*y2)+15*x2*x2);
# 			}break;
# 			default:
# 			err = 1;
# 		}
# 		*dx *= norm[n];
# 		*dy *= norm[n];
# 		return(err);
# 	}

    return d
