import numpy as np
def makeMask( m, masklist, stampD, maskScalingFactor ):

    m.pMask = np.ones( stampD, dtype=int )
    m.cMask = m.pMask

    xSensor, ySensor = np.mgrid[-(stampD/2-0.5):(stampD/2+0.5),-(stampD/2-0.5):(stampD/2+0.5)]

    #print xSensor.shape, ySensor.shape

    xSensor=xSensor/(stampD/2/m.sensorFactor)
    ySensor=ySensor/(stampD/2/m.sensorFactor)
    rMask = m.apertureDiameter / (2 * m.focalLength/m.offset)*maskScalingFactor

    #print "rMask",rMask
    #print masklist.shape[0]
#    print masklist
    for ii in range(masklist.shape[0]):
#        print "II",ii

        r = np.sqrt( (xSensor-masklist[ii,0])**2 + (ySensor-masklist[ii,1])**2 )

        # Initialize both mask elements to the opposite of the pass/block boolean
        pMask = (1-masklist[ii,3])*np.ones( (stampD,stampD), dtype=int )
        cMask = (1-masklist[ii,3])*np.ones( (stampD,stampD) )
#        print masklist[ii,3]
#        print pMask
#        print pMask.shape
        # Find the indices that correspond to the mask element, set them to the pass/block boolean
        idx =  r <= masklist[ii,2]
        if (masklist[ii,3] >= 1):
            aidx = np.nonzero(r <= masklist[ii,2]*(1+m.boundaryT*m.pixelSize/rMask)) #Bo: make a mask >r so that we can keep a larger area of S
        else:
            aidx = np.nonzero(r <= masklist[ii,2]*(1-m.boundaryT*m.pixelSize/rMask))
        pMask[idx] = masklist[ii,3]
        cMask[aidx] = masklist[ii,3]

        # Multiplicatively add the current mask elements to the model masks
        m.pMask = m.pMask*pMask    #padded mask - for use at the offset planes
        m.cMask = m.cMask*cMask    #non-padded mask corresponding to aperture

        pMask = m.pMask
        cMask = m.cMask

    #print pMask.max()
    #print cMask.max()
    return pMask, cMask

