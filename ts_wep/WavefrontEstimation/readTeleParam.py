import numpy as np
import os
import string
from interpMaskParam import interpMaskParam


def readTeleParam(m,instruFile,I1fldx,I1fldy,I2fldx,I2fldy):

    I1fldr = np.sqrt(I1fldx**2+I1fldy**2)
    I2fldr = np.sqrt(I2fldx**2+I2fldy**2)
    fid=open(os.path.join('',instruFile))

    iscomment = False
    for line in fid:
        line = line.strip()
        if (line.startswith('###')):
            iscomment = ~iscomment
        if ((line.startswith('#') == False) and (iscomment == False) and len(line) >0):
            if (line.startswith('Obscuration')):
                m.obscuration = string.atof(line.split()[-1])
            if (line.startswith('Focal_length')):
                m.focalLength = string.atof(line.split()[-1])
            if (line.startswith('Aperture_diameter')):
                m.apertureDiameter = string.atof(line.split()[-1])
            if (line.startswith('Offset')):
                m.offset = string.atof(line.split()[-1])
            if (line.startswith('Pixel_size')):
                m.pixelSize = string.atof(line.split()[-1])
            if (line.startswith('Mask_param')):
                maskParam = os.path.join('../conf/',line.split()[-1])
            if (line.startswith('Marginal_fl')):
                m.marginalFL = string.atof(line.split()[-1])

    if (not (I1fldx==0 and I1fldy==0 )):
        m.I1maskCa, m.I1maskRa, m.I1maskCb, m.I1maskRb = interpMaskParam(maskParam,I1fldr)

    if ( not (I2fldx==0 and I2fldy==0 )):
        m.I2maskCa, m.I2maskRa, m.I2maskCb, m.I2maskRb = interpMaskParam(maskParam,I2fldr)

    fid.close()                
    
    return m



