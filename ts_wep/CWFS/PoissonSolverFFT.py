

from createSignal import createSignal
from ZernikeMaskedFit import ZernikeMaskedFit
from extractArray import extractArray
from padArray import padArray

import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage

def PoissonSolverFFT(m, cliplevel):
    '''Poisson Solver using an FFT

    inputs
    outputs
    '''
    aperturePixelSize = (m.apertureDiameter*m.sensorFactor/m.sensorSamples)
    # print "FFT"
    #    kk = np.arange(-0.5,0.5-1/m.padDim,1/m.padDim) / aperturePixelSize
    v, u = np.mgrid[-0.5/aperturePixelSize:(0.5)/aperturePixelSize:1/m.padDim/aperturePixelSize,-0.5/aperturePixelSize:(0.5)/aperturePixelSize:1/m.padDim/aperturePixelSize]

    u2v2 = -4*(np.pi**2)*(u*u + v*v)

    # Set origin to Inf and 0 to result in 0 at origin after filtering
    ctrIdx = np.floor(m.padDim/2)
    u2v2[ctrIdx, ctrIdx] = np.inf

    m = createSignal( m , cliplevel)

    s=extractArray(m.S,m.sensorSamples) #test

#    fig = plt.figure(figsize=(5, 2.2))
#    fig.subplots_adjust(left=0.12, right=0.95, bottom=0.2, top=0.9,
#                        hspace=0.01, wspace=0.01)
#    ax = plt.subplot(111, aspect='equal')
#    ax.imshow(s.T, origin='lower')
#    plt.show()

    #find the indices for a ring of pixels just ouside and just inside the
    #aperture for use in setting dWdn = 0

    struct = ndimage.generate_binary_structure(2, 1)
    struct = ndimage.iterate_structure(struct, m.boundaryT)
    #print struct
    ApringOut = np.logical_xor(ndimage.morphology.binary_dilation(m.pMask,structure=struct),m.pMask).astype(int)
    ApringIn = np.logical_xor(ndimage.morphology.binary_erosion(m.pMask,structure=struct),m.pMask).astype(int)
#    print "ApringOut",ApringOut,np.nonzero(ApringOut)#out is different in is the same
    bordery, borderx = np.nonzero(ApringOut) # x= y and y=x

    if (m.compMode == 'zer'):
        zc = np.zeros((m.numTerms,m.innerItr))
#        print "ZC ONE",zc.shape

    #*************************************************************************
    #initial BOX 3 - put signal in boundary (since there's no existing
    #Sestimate, S just equals m.S
    S = m.S.copy()

    for jj in range(int(m.innerItr)):

        #*************************************************************************
        #BOX 4 - forward filter: forward FFT, divide by u2v2, inverse FFT
        SFFT = np.fft.fftshift(np.fft.fft2( np.fft.fftshift(S) ))
        #print SFFT.shape, u2v2.shape
        W = np.fft.fftshift( np.fft.irfft2( np.fft.fftshift( SFFT/u2v2 ), s=S.shape) )

        #*************************************************************************
        #BOX 5 - Wavefront estimate (includes zeroing offset & masking to the aperture size)
        West = extractArray( W, m.sensorSamples )
        #print "WEST", West.shape, W.shape

        offset = West[m.pMask==1].mean()
        West = West - offset
        West[m.pMask == 0] = 0

        if (m.compMode =='zer'):
            ySensor, xSensor = np.mgrid[-(m.sensorSamples/2-0.5):(m.sensorSamples/2+0.5),-(m.sensorSamples/2-0.5):(m.sensorSamples/2+0.5)]
            xSensor=xSensor/(m.sensorSamples/2/m.sensorFactor)
            ySensor=ySensor/(m.sensorSamples/2/m.sensorFactor)
#            print "ZC",zc.shape
            zc[:,jj] = ZernikeMaskedFit(West, xSensor, ySensor, m.numTerms, m.pMask, m.zobsR )
#            print "ZC O",zc.shape

        #*************************************************************************
        #BOX 6 - set dWestimate/dn = 0 around boundary
        WestdWdn0 = West.copy()

        #do a 3x3 average around each border pixel, including only those pixels inside the aperture
#        print borderx
#        print borderx.shape
        for ii in range(len(borderx)):
            reg = West[borderx[ii]-m.boundaryT:borderx[ii]+m.boundaryT+1,bordery[ii]-m.boundaryT:bordery[ii]+m.boundaryT+1]
            intersectIdx = ApringIn[borderx[ii]-m.boundaryT:borderx[ii]+m.boundaryT+1,bordery[ii]-m.boundaryT:bordery[ii]+m.boundaryT+1]
            WestdWdn0[borderx[ii],bordery[ii]] = reg[np.nonzero(intersectIdx)].mean()

        #*************************************************************************
        #BOX 7 - Take Laplacian to find sensor signal estimate

        Wxx = np.zeros((m.sensorSamples,m.sensorSamples))
        Wyy = np.zeros((m.sensorSamples,m.sensorSamples))
        Wt = WestdWdn0.copy()
#        print Wxx.shape, Wt.shape
#        print Wxx[:,1:-1].shape
        Wxx[:,1:-1] = (Wt[:,0:-2]-2*Wt[:,1:-1]+Wt[:,2:])/aperturePixelSize**2
        Wyy[1:-1,:] = (Wt[0:-2,:]-2*Wt[1:-1,:]+Wt[2:,:])/aperturePixelSize**2
        del2W = Wxx + Wyy
        Sest = padArray(del2W, m.padDim)

        #*************************************************************************
        #BOX 3 - Put signal back inside boundary, leaving the rest of Sestimate
        Sest[m.pMaskPad==1] = m.S[m.pMaskPad==1]
        S = Sest


    m.West = West.copy()
    if (m.compMode == 'zer'):
        m.zc = zc

    return m
