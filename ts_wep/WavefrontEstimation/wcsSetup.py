from readTeleParam import readTeleParam
from readAlgoParam import readAlgoParam
from upResolution import upResolution
from makeCornerMasklist import makeCornerMasklist
from makeMask import makeMask
from padArray import padArray

import numpy as np

class DataStruct(object):
    pass
    
def wcsSetup(I1,I1fldx,I1fldy,I2,I2fldx,I2fldy,instruFile,algoFile):
    m = DataStruct()
    m.I1=I1;
    m.I2=I2;
    m = readTeleParam(m,instruFile,I1fldx,I1fldy,I2fldx,I2fldy);
    m = readAlgoParam(m,algoFile,I1fldx,I1fldy,I2fldx,I2fldy);

    m.fno=m.focalLength/m.apertureDiameter
    m.sensorSamples = I1.shape[0]
#    print m.sensorSamples
    if 'PoissonSolver' in m.__dict__.keys():
#        print m.PoissonSolver
        if ((m.PoissonSolver == 'fft') and (m.padDim == 999)):
            m.padDim = 2**np.ceil(np.log2(m.sensorSamples))
#            print m.padDim

    m.sensorFactor = m.sensorSamples / (m.offset * m.apertureDiameter / m.focalLength/ m.pixelSize)
    m.sensorWidth = (m.apertureDiameter*m.offset/m.focalLength) * m.sensorFactor
    m.donutR = m.pixelSize*(m.sensorSamples/m.sensorFactor)/2
    if 'marginalFL' not in m.__dict__.keys():
        m.maskScalingFactor = 1
    else:
        m.maskScalingFactor = m.focalLength/m.marginalFL

    if (('upReso'  in m.__dict__.keys()) and ( m.upReso > 1)):
        # TODO not tested
        m.I1 = upResolution(m.I1,m.upReso,m.sensorSamples*m.upReso,m.sensorSamples*m.upReso);
        m.I2 = upResolution(m.I2,m.upReso,m.sensorSamples*m.upReso,m.sensorSamples*m.upReso);
        m.pixel_size=m.pixel_size/m.upReso;
        m.sensorSamples = m.sensorSamples*m.upReso;

    if (I1fldx==0 and I1fldy==0 ):
        if m.obscuration == 0:
            m.I1masklist = np.array([0, 0, 1, 1])
        else:
            m.I1masklist = np.array([[0, 0, 1, 1], [0, 0, m.obscuration, 0]])
    else:
        m.I1masklist=makeCornerMasklist(m.obscuration,m.I1maskCa,m.I1maskRa,m.I1maskCb,m.I1maskRb,I1fldx,I1fldy)

    if (I2fldx==0 and I2fldy==0 ):
        if m.obscuration == 0:
            m.I2masklist = np.array([0, 0, 1, 1])
        else:
            m.I2masklist = np.array([[0, 0, 1, 1], [0, 0, m.obscuration, 0]])
    else:
        m.I2masklist=makeCornerMasklist(m.obscuration,m.I2maskCa,m.I2maskRa,m.I2maskCb,m.I2maskRb,I2fldx,I2fldy)

#    print  m.I2masklist
    m.I1pMask, m.I1cMask = makeMask(m,m.I1masklist,m.sensorSamples,1)
    m.I2pMask, m.I2cMask = makeMask(m,m.I2masklist,m.sensorSamples,1)
    m.pMask = m.I1pMask*m.I2pMask 
    m.cMask = m.I1cMask*m.I2cMask


    if (('PoissonSolver'  in m.__dict__.keys()) and (m.PoissonSolver == 'fft')):
        m.pMaskPad = padArray( m.pMask, m.padDim );
        m.cMaskPad = padArray( m.cMask, m.padDim );

    return m
