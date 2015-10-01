import numpy as np
import os
import string
from interpMaskParam import interpMaskParam
from getOffAxisCorr import getOffAxisCorr
from numpy import loadtxt

def readAlgoParam(m,algoFile,I1fldx,I1fldy,I2fldx,I2fldy):

    fid=open(os.path.join('',algoFile))

    iscomment = False
    for line in fid:
        line = line.strip()
        if (line.startswith('###')):
            iscomment = ~iscomment
        if ((line.startswith('#') == False) and (iscomment == False) and len(line) >0):
            if (line.startswith('PoissonSolver')):
                m.PoissonSolver = line.split()[1]
            if (line.startswith('Num_of_Zernikes')):
                m.numTerms = string.atoi(line.split()[1])
            if (line.startswith('ZTerms')):
                m.ZTerms = map(int,line.split()[1:])
            if (line.startswith('Num_of_outer_itr')):
                m.outerItr = string.atoi(line.split()[1])
            if (line.startswith('Num_of_inner_itr')):
                m.innerItr = string.atoi(line.split()[1])
            if (line.startswith('Zernikes')):
                m.zobsR = string.atoi(line.split()[1])
            if (line.startswith('Increase_resolution')):
                m.upReso = string.atof(line.split()[1])
            if (line.startswith('FFT_dimension')):
                m.padDim = string.atof(line.split()[2])
            if (line.startswith('Feedback_gain')):
                m.feedbackGain = string.atof(line.split()[1])
            if (line.startswith('Compensator_oversample')):
                m.compOversample = string.atof(line.split()[1])
            if (line.startswith('Compensator_mode')):
                m.compMode = line.split()[1]
            if (line.startswith('OffAxis_poly_order')):
                m.offAxisPolyOrder = string.atoi(line.split()[1])
            if (line.startswith('Boundary_thickness')):
                m.boundaryT = string.atoi(line.split()[2])
            if (line.startswith('Compensation_sequence')):
                m.compSequence = loadtxt(os.path.join('../conf/',line.split()[1]))
            if (line.startswith('Sumclip_sequence')):
                m.sumclipSequence = loadtxt(os.path.join('../conf/',line.split()[1]))
            if (line.startswith('Image_formation')):
                m.imageFormation = line.split()[1]
            if (line.startswith('Minimization')):
                m.minimization = line.split()[1]
    fid.close()

    if  not ('ZTerms' in m.__dict__.keys()):
        m.ZTerms=np.arange(m.numTerms)+1 #starts from 1

    if 'zobsR' in m.__dict__.keys():
        if ((m.zobsR==1) and ('obscuration' in m.__dict__.keys())):
            m.zobsR=m.obscuration

    #if m.outerItr is large, and m.compSequence is too small, we fill in the
    #rest in m.compSequence
    #print m.compSequence.shape[0]
    if (m.compSequence.shape[0] < m.outerItr): 
        if (len(m.compSequence.shape) == 1):
            #TODO resize compSequence to be m.outer and set all etra values to compSequence[-1]
            m.compSequence[m.compSequence.shape[0]+1:m.outerItr] = m.compSequence[-1]
        else:
            #TODO for all dimensions resize compSequence to be m.outer and set all etra values to 1
            m.compSequence[:,m.compSequence.shape[1]+1:m.outerItr] = 1

    I1fldr = np.sqrt(I1fldx**2+I1fldy**2)
    I2fldr = np.sqrt(I2fldx**2+I2fldy**2)

    if 'offAxisPolyOrder' in m.__dict__.keys():
        m.I1offAxis_coeff = np.zeros((4,(m.offAxisPolyOrder+1)*(m.offAxisPolyOrder+2)/2))
        m.I1offAxis_coeff[0,:] = getOffAxisCorr('../conf/lsst_offAxis_cxin_poly%d.txt'%(m.offAxisPolyOrder),I1fldr)
        m.I1offAxis_coeff[1,:] = getOffAxisCorr('../conf/lsst_offAxis_cyin_poly%d.txt'%(m.offAxisPolyOrder),I1fldr)
        m.I1offAxis_coeff[2,:] = getOffAxisCorr('../conf/lsst_offAxis_cxex_poly%d.txt'%(m.offAxisPolyOrder),I1fldr)
        m.I1offAxis_coeff[3,:] = getOffAxisCorr('../conf/lsst_offAxis_cyex_poly%d.txt'%(m.offAxisPolyOrder),I1fldr)

        m.I2offAxis_coeff = np.zeros((4,(m.offAxisPolyOrder+1)*(m.offAxisPolyOrder+2)/2))
        m.I2offAxis_coeff[0,:] = getOffAxisCorr('../conf/lsst_offAxis_cxin_poly%d.txt'%(m.offAxisPolyOrder),I2fldr)
        m.I2offAxis_coeff[1,:] = getOffAxisCorr('../conf/lsst_offAxis_cyin_poly%d.txt'%(m.offAxisPolyOrder),I2fldr)
        m.I2offAxis_coeff[2,:] = getOffAxisCorr('../conf/lsst_offAxis_cxex_poly%d.txt'%(m.offAxisPolyOrder),I2fldr)
        m.I2offAxis_coeff[3,:] = getOffAxisCorr('../conf/lsst_offAxis_cyex_poly%d.txt'%(m.offAxisPolyOrder),I2fldr)

    return m
