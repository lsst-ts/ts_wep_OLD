from astropy.io import fits
from imageCoCenter import imageCoCenter
from PoissonSolverFFT import PoissonSolverFFT
from PoissonSolverExp import PoissonSolverExp
from compensate import compensate
from ZernikeEval import ZernikeEval
from ZernikeAnnularEval import ZernikeAnnularEval
import copy

import numpy as np
from wcsSetup import wcsSetup

def wcs(I1File, I1fldx, I1fldy, I2File, I2fldx, I2fldy, instruFile, algoFile, model):

    if (isinstance(I1File, str) and isinstance(I2File,str)):
        if (I1File.endswith(".txt") or I2File.endswith(".TXT")):
            I1 = np.loadtxt(I1File) 	    
            I2 = np.loadtxt(I2File)
            # flip along horizontal axis #np.filpud() also works
            # because of how ZEMAX writes out#images and how MATLAB reads in images
            I1 = I1[ ::-1,:] 
            I2 = I2[ ::-1,:] 
        else:
            I1HDU = fits.open(I1File)
            I1 = I1HDU[0].data
            I2HDU = fits.open(I2File)
            I2 = I2HDU[0].data
            I1HDU.close()
            I2HDU.close()
    else:
        I1 = I1File    
        I2 = I2File    

# MATLAB Image convolution (for investigation of pixel aliasing effect)
# convkernel=load('test/gau_6x6_fwhm1.txt');
# convkernel=load('test/gau_10x10_fwhm3.txt');
# convkernel=load('test/gau_40x40_fwhm10.txt');
# convkernel=load('test/gau_100x100_fwhm30.txt');
# I1= conv2(I1, convkernel, 'same');
# I2= conv2(I2, convkernel, 'same');
# I1=downResolution(I1,10,120,120);
# I2=downResolution(I2,10,120,120);
    m = wcsSetup(I1, I1fldx, I1fldy, I2, I2fldx, I2fldy, instruFile, algoFile)

    m.converge = np.zeros((m.numTerms, m.outerItr + 1))

    #for estimating m.Wconverge
    ySensor, xSensor = np.mgrid[-(m.sensorSamples/2-0.5):(m.sensorSamples/2+0.5), \
                                -(m.sensorSamples/2-0.5):(m.sensorSamples/2+0.5)]
    xSensor=xSensor/(m.sensorSamples/2/m.sensorFactor)
    ySensor=ySensor/(m.sensorSamples/2/m.sensorFactor)
    r2Sensor=xSensor**2+ySensor**2
    idx=(r2Sensor>1) | (r2Sensor<m.obscuration**2)
    xSensor[idx]=np.nan
    ySensor[idx]=np.nan

    m = imageCoCenter(m, I1fldx, I1fldy, I2fldx, I2fldy)        
    #print m.__dict__.keys()

    m0 = copy.deepcopy(m)
    
    if m.compMode == 'zer':
        ztot = np.zeros(m.numTerms)    

        m.caustic = 0
        if 'Axis' in model:        #onAxis or offAxis
            m.I1, tmp = compensate(m0, ztot, 'intra', 1, I1fldx, I1fldy, model)        
            m.I2, tmp = compensate(m0, ztot, 'extra', 1, I2fldx, I2fldy, model)        

        #     model='paraxial';
        # matrix or element multiplication
        if (I1fldx != I2fldx or I1fldy != I2fldy):        
            m.I1 = m.I1 * m.pMask            
            m.I2 = m.I2 * np.rot90(m.pMask, 2)
            m.I1 = m.I1 / np.sum(m.I1)            
            m.I2 = m.I2 / np.sum(m.I2)        #no need vignetting correction, this is after masking already

        if m.PoissonSolver == 'fft':
            m = PoissonSolverFFT(m, m.sumclipSequence[0])
#            print m.converge.shape, m.zc.shape, m.innerItr
            m.converge[:, 0] = ztot + m.zc[:, m.innerItr-1]
        else:        
            m = PoissonSolverExp(m)        
            m.converge[:, 0] = ztot + m.zc        

        #    m.West includes Zernikes presented by m.zc
        m.Wconverge=m.West;
                    
        for j in range(int(m.outerItr)):
            if not m.caustic:            
                if (m.PoissonSolver == 'fft'):
                    ztmp = m.zc[:, -1]
                else:                
                    ztmp = m.zc                
#                print m.compSequence.shape
                if (m.compSequence.ndim == 1):                
                    ztmp[m.compSequence[j]:] = 0                
                else:                
                    ztmp = ztmp * m.compSequence[:, j]

                ztot = ztot + ztmp * m.feedbackGain            
                
                m.I1, caustic = compensate(m0, ztot,'intra', 1, I1fldx, I1fldy, model)            
                if caustic > 0:                    
                    m.caustic = caustic

                m.I2, caustic = compensate(m0, ztot, 'extra', 1, I2fldx, I2fldy, model)            
                if caustic > 0:                    
                    m.caustic = caustic
                if (I1fldx != I2fldx or I1fldy != I2fldy):                
                    m.I1 = m.I1 * m.pMask
                    m.I2 = m.I2 * np.rot90(m.pMask, k=2);                
                    m.I1 = m.I1 / np.sum(m.I1)                    
                    m.I2 = m.I2 / np.sum(m.I2)  #no need vignetting correction, this is after masking already

                if (m.PoissonSolver == 'fft'):                
                    m = PoissonSolverFFT(m, m.sumclipSequence[j])                
                    m.converge[:, j+1] = ztot + m.zc[:, -1] 
                else:                
                    m = PoissonSolverExp(m)                
                    m.converge[:, j+1] = ztot + m.zc
                    
                #m.West is the estimated wavefront from the last run of Poisson
                #solver. ztot is what had be compensated before that run.
                #m.West includes two parts: latest m.zc, and m.Wres
                #m.Wres is the residual wavefront on top of m.converge(:,end),
                #m.Wres is only available for the iterative FFT algorithm.
                if (m.zobsR==0):
                    m.Wconverge=ZernikeEval(np.concatenate(([0,0,0],ztot[3:]),axis=1),\
                                            xSensor,ySensor)+m.West
                else:
                    m.Wconverge=ZernikeAnnularEval(np.concatenate(([0,0,0],ztot[3:]),axis=1),\
                                                   xSensor,ySensor,m.zobsR)+m.West;
                    
            else:            # once we run into caustic, stop here, results may be
                             # close to real aberration. Continuation may lead to disatrous results
                m.converge[:, j+1] = m.converge[:, j]

    elif (m.compMode == 'opd'):

        wtot = np.zeros(m.sensorSamples, m.sensorSamples)    
        
        m.caustic = 0    
        if ('Axis' in model):        #onAxis or offAxis
            compensate(m0, wtot, 'intra', 1, I1fldx, I1fldy, model)        
            compensate(m0, wtot, 'extra', 1, I2fldx, I2fldy, model)        

        if (I1fldx != I2fldx or I1fldy != I2fldy):        
            m.I1 = m.I1 * m.pMask            
            m.I2 = m.I2 * np.rot90(m.pMask, k=2)        
            m.I1 = m.I1 / np.sum(m.I1)            
            m.I2 = m.I2 / np.sum(m.I2)    #no need vignetting correction, this is after masking already

        if (m.PoissonSolver == 'fft'):
            m = PoissonSolverFFT(m, m.sumclipSequence(1))        
        else:        
            m = PoissonSolverExp(m)        

        Wconverge = wtot + m.West    
        m.converge[:, 0] = ZernikeMaskedFit(Wconverge, xSensor, ySensor, m.numTerms, m.pMask, m.zobsR)    
        for j in range(int(m.outerItr)):
            if not m.caustic: 
                wtmp = m.West            
                wtot = wtot + wtmp * m.feedbackGain            
                
                m.I1, caustic = compensate(m0, wtot, 'intra', 1, I1fldx, I1fldy, model)            
                if caustic > 0:                    
                    m.caustic = caustic

                m.I2, caustic = compensate(m0, wtot, 'extra', 1, I2fldx, I2fldy, model)            
                if caustic > 0:                    
                    m.caustic = caustic

                if (I1fldx != I2fldx or I1fldy != I2fldy):                
                    m.I1 = m.I1 * m.pMask                    
                    m.I2 = m.I2 * np.rot90(m.pMask, k=2) 
                    m.I1 = m.I1 / np.sum(m.I1)                    
                    m.I2 = m.I2 / np.sum(m.I2)   #no need vignetting correction this is after masking already

                if (m.PoissonSolver == 'fft'):
                    m = PoissonSolverFFT(m, m.sumclipSequence[j])                
                else:                
                    m = PoissonSolverExp(m)                

                Wconverge = wtot + m.West            
                m.converge[:, j] = ZernikeMaskedFit(Wconverge, xSensor, ySensor, m.numTerms,
                                                    m.pMask, m.zobsR)            
            else:         #once we run into caustic, stop here, results may be close to real aberration. Continuation may lead to disatrous results
                m.converge[:, j] = m.converge[:, j]            
    m.I1 = I1
    m.I2 = I2
    
    return m


