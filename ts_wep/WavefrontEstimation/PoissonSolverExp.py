from getdIandI import getdIandI
from ZernikeEval import ZernikeEval
from ZernikeGrad import ZernikeGrad
from ZernikeAnnularEval import ZernikeAnnularEval
from ZernikeAnnularGrad import ZernikeAnnularGrad
import numpy as np

def PoissonSolverExp(m):

    m = getdIandI(m)
    ySensor, xSensor = np.mgrid[-(m.sensorSamples/2-0.5):(m.sensorSamples/2+0.5), \
                                -(m.sensorSamples/2-0.5):(m.sensorSamples/2+0.5)]
    
    xSensor=xSensor/(m.sensorSamples/2/m.sensorFactor)*m.cMask
    ySensor=ySensor/(m.sensorSamples/2/m.sensorFactor)*m.cMask

    F=np.zeros(m.numTerms)
    dZidx=np.zeros((m.numTerms,m.sensorSamples,m.sensorSamples))
    dZidy=dZidx.copy()

    aperturePixelSize = (m.apertureDiameter*m.sensorFactor/m.sensorSamples)
    zcCol = np.zeros(m.numTerms)
    for i in range(int(m.numTerms)):
        zcCol[i]=1;
        #we integrate, instead of decompose, integration is faster. Also, decomposition is ill-defined on m.cMask. 
        #Using m.pMask, the two should give same results.
        if (m.zobsR>0):
            F[i]=np.sum(np.sum(m.dI*ZernikeAnnularEval(zcCol,xSensor,ySensor,
                                                        m.zobsR)))*aperturePixelSize**2
            dZidx[i,:,:] = ZernikeAnnularGrad(zcCol,xSensor,ySensor,m.zobsR,'dx')
            dZidy[i,:,:] = ZernikeAnnularGrad(zcCol,xSensor,ySensor,m.zobsR,'dy')
        else:
            F[i]=np.sum(np.sum(m.dI*ZernikeEval(zcCol,xSensor,ySensor)))*aperturePixelSize**2
            dZidx[i,:,:]=ZernikeGrad(zcCol,xSensor,ySensor,'dx')
            dZidy[i,:,:]=ZernikeGrad(zcCol,xSensor,ySensor,'dy')
        zcCol[i]=0

    m.Mij=np.zeros((m.numTerms,m.numTerms))
    for i in range(m.numTerms):
        for j in range(m.numTerms):
            m.Mij[i,j]=aperturePixelSize**2/(m.apertureDiameter/2)**2* \
                        np.sum(np.sum(m.I*(dZidx[i,:,:].squeeze()*dZidx[j,:,:].squeeze() \
                                           + dZidy[i,:,:].squeeze()*dZidy[j,:,:].squeeze() )))

    dz=2*m.focalLength*(m.focalLength-m.offset)/m.offset
    m.zc = np.zeros(m.numTerms)
    idx = [x-1 for x in m.ZTerms]
    zc_tmp = np.dot(np.linalg.pinv(m.Mij[:,idx][idx]),F[idx])/dz #phi in GN paper is phase, phi/(2pi)*lambda=W
    m.zc[idx]=zc_tmp

    if (m.zobsR>0):
        m.West=ZernikeAnnularEval(np.concatenate(([0,0,0],m.zc[3:]),axis=1),xSensor,ySensor,m.zobsR)
    else:
        m.West=ZernikeEval(np.concatenate(([0,0,0],m.zc[3:]),axis=1),xSensor,ySensor)

    return m
