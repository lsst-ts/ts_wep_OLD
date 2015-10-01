import numpy as np
import scipy.interpolate as interpolate
from aperture2image import aperture2image
from showProjection import showProjection
from getCenterAndR_ef import getCenterAndR_ef
from padArray import padArray
from extractArray import extractArray
import scipy.ndimage as ndimage

def compensate(m, zcCol, atype, oversample,fldx,fldy,model):

    #print zcCol.ndim, zcCol.shape[0], m.numTerms
    if ((zcCol.ndim ==1) and (len(zcCol) !=m.numTerms)):
        raise Exception('input:size','zcCol in compensate needs to be a %d row column vector\n'%m.numTerms)

    if (('nc' in model) and ('offAxis' in model)):
        if (zcCol.ndim ==1):
            phaseScreen = np.loadtxt('wcs/conf/lsst_NE_screen_annuZ.txt')
        else:
             phaseScreen = np.loadtxt('wcs/conf/lsst_NE_screen_opd.txt')
        zcCol = zcCol + m.alambda * phaseScreen

    if (atype == 'intra'):
        myIp = m.I1
    elif (atype == 'extra'):
        myIp = m.I2

    sm, sn = myIp.shape
    if (sm !=sn ):
        raise Exception('=========Error in fcompensate.m: real image is not a square. \n')
        return
    
    projSamples = sm*oversample

    # Let us create a look-up table for x -> xp first. 
    luty, lutx = np.mgrid[-(projSamples/2-0.5):(projSamples/2+0.5),-(projSamples/2-0.5):(projSamples/2+0.5)]
    lutx = lutx/(projSamples/2/m.sensorFactor)
    luty= luty/(projSamples/2/m.sensorFactor)

    #set up the mapping
    # [lutxp lutyp J]=aperture2image_pMask(m,zcCol,lutx,luty,m.pMask,atype,fldx,fldy,model);
    lutxp, lutyp, J = aperture2image(m,zcCol,lutx,luty,projSamples,atype,fldx,fldy,model)
#    print "J",J.shape

    show_lutxyp = showProjection(lutxp, lutyp, m.sensorFactor, projSamples,0);

    realcx, realcy, tmp = getCenterAndR_ef(myIp)
    show_lutxyp=padArray(show_lutxyp,projSamples+20)

    struct0 = ndimage.generate_binary_structure(2, 1)
    struct = ndimage.iterate_structure(struct0, 4)
    struct=ndimage.morphology.binary_dilation(struct,structure=struct0)
    struct=ndimage.morphology.binary_dilation(struct,structure=struct0).astype(int)
    show_lutxyp = ndimage.morphology.binary_dilation(show_lutxyp,structure=struct)
    show_lutxyp = ndimage.morphology.binary_erosion(show_lutxyp,structure=struct)
    show_lutxyp=extractArray(show_lutxyp,projSamples)

    projcx, projcy, tmp = getCenterAndR_ef(show_lutxyp.astype(float))
    projcx = projcx/(oversample)
    projcy = projcy/(oversample)

    
    # +(-) means we need to move image to the right (left)
    shiftx = (projcx - realcx)
    # +(-) means we need to move image upward (downward)
    shifty = (projcy - realcy)
    myIp= np.roll(myIp, int(np.round(shifty)), axis=0)
    myIp= np.roll(myIp, int(np.round(shiftx)), axis=1)

    #% let's construct the interpolant, to get the intensity on (x',p') plane
    #% that corresponds to the grid points on (x,y)
#    print np.mgrid[-((sm)/2-0.5):(sm/2+0.5)]
    yp,xp = np.mgrid[-(sm/2-0.5):(sm/2+0.5),-(sm/2-0.5):(sm/2+0.5)]
    xp = xp/(sm/2/m.sensorFactor)
    yp= yp/(sm/2/m.sensorFactor)

    # xp = reshape(xp,sm^2,1);
    # yp = reshape(yp,sm^2,1);
    # myIp = reshape(myIp,sm^2,1);
    #  
    # FIp = TriScatteredInterp(xp,yp,myIp,'nearest');
    # lutIp = FIp(lutxp, lutyp);
 
    lutxp[np.isnan(lutxp)] = 0
    lutyp[np.isnan(lutyp)] = 0
 
#    lutIp=interp2(xp,yp,myIp,lutxp,lutyp)
#    print xp.shape, yp.shape, myIp.shape
#    print lutxp.ravel()
#    print xp[:,0],yp[0,:]
    ip = interpolate.RectBivariateSpline( yp[:,0], xp[0,:], myIp,kx=1,ky=1)

#    ip = interpolate.interp2d(xp, yp, myIp)
#    ip = interpolate.interp2d(xp, yp, myIp)
#    print lutxp.shape, lutyp.shape
#    lutIp = ip(0.5, -0.5)
#    print lutIp, 'lutIp1'
#    lutIp = ip([-0.1],[-0.1])
#    print lutIp, 'lutIp2'
#    lutIp = ip(np.array(0.5,-0.1), np.array(-0.5, -0.1))
#    print lutIp, 'lutIp12',lutxp.ravel()[0:10]
    lutIp = np.zeros(lutxp.shape[0]*lutxp.shape[1])
    for i,(xx,yy) in enumerate(zip(lutxp.ravel(),lutyp.ravel())):
        lutIp[i] = ip(yy, xx)
    lutIp = lutIp.reshape(lutxp.shape)
    
    myI = lutIp * J

    

    if (atype == 'extra'):
        myI = np.rot90( myI, k=2 )


    #if we want the compensator to drive down tip-tilt
    # myI = offsetImg(-shiftx, -shifty, myI);
    # myI=circshift(myI,[round(-shifty) round(-shiftx)]);

    myI[np.isnan(myI)]=0
    # myI < 0 will not be physical, remove that region
    # x(myI<0) = NaN;
    caustic = 0;
    if (myI.min() < 0):
        print('WARNING: negative scale parameter, image is within caustic, zcCol (in um)=\n');

#    for i in range(len(zcCol)):
#        print zcCol[i]
#        print('%5.2f '%(zcCol[i]*1.e6))
#    print('\n');
        caustic=1;

    myI[myI<0] = 0;
    if (oversample>1):
        newI=downResolution(myI,oversample,sm,sn)
    else:
        newI=myI

    return newI, caustic
