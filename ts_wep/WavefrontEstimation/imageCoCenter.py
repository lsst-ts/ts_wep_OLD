from getCenterAndR_ef import getCenterAndR_ef
import numpy as np

def imageCoCenter(m,I1fldx,I1fldy,I2fldx,I2fldy):

    x1, y1, tmp = getCenterAndR_ef(m.I1)
    x2, y2, tmp = getCenterAndR_ef(m.I2)
#    print x1,y1,x2,y2,tmp
    stampCenterx1=m.sensorSamples/2. +0.5
    stampCenterx2=m.sensorSamples/2. +0.5   
    stampCentery1=m.sensorSamples/2. +0.5
    stampCentery2=m.sensorSamples/2. +0.5 
    radialShift=3.5*m.upReso*(m.offset/1e-3)*(10e-6/m.pixelSize)
    
    I1fldr = np.sqrt(I1fldx**2+I1fldy**2);
    I2fldr = np.sqrt(I2fldx**2+I2fldy**2);
    I1radialShift = radialShift*I1fldr/1.75
    I2radialShift = radialShift*I2fldr/1.75
    if (I1fldr>1.75): 
        I1radialShift=0

    if (I2fldr>1.75):
        I2radialShift=0

    if I1fldr != 0:
        I1c = I1fldx/I1fldr
        I1s = I1fldy/I1fldr
    else:
        I1c=0
        I1s=0

    if I2fldr !=0:
        I2c = I2fldx/I2fldr
        I2s = I2fldy/I2fldr
    else:
        I2c=0
        I2s=0
    
    stampCenterx1=stampCenterx1+I1radialShift*I1c
    stampCentery1=stampCentery1+I1radialShift*I1s
    stampCenterx2=stampCenterx2-I2radialShift*I2c
    stampCentery2=stampCentery2-I2radialShift*I2s

 #   print "INDEX",y1 #roll is correct 
    m.I1=np.roll(m.I1,int(np.round(stampCentery1-y1)),axis=0)
    m.I1=np.roll(m.I1,int(np.round(stampCenterx1-x1)),axis=1)
    m.I2=np.roll(m.I2,int(np.round(stampCentery2-y2)),axis=0)
    m.I2=np.roll(m.I2,int(np.round(stampCenterx2-x2)),axis=1)

    
    return m
