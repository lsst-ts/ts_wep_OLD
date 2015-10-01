import numpy as np
import scipy as sci
from padArray import padArray

def opd2psf(opd,imagedelta,sensorFactor,fno,wavelength):
    """OPD 
    wavefront OPD in wave
    imagedelta in micron
    wavelength in meter
    """
    m, n = opd.shape

    if (m != n):
        print 'warning: opd is not a square stamp'
        return 0

    k= 1/imagedelta*fno*wavelength/1e-6
    padding=k/sensorFactor
    if (padding < 1):
        print 'opd2psf: sampling too low, data inaccurate'
        return 0

    sensorSamples = np.max(m,n)
    pupil = opd.copy().tocsr()
    pupil.data.fill(1).todense()

    opd[np.isnan(opd)] = 0
    N = np.round(padding*sensorSamples)
    pupil = padArray(pupil,N)

    opd = padArray(opd,N)

    z = pupil * np.exp(-2*pi*1j*opd)
    z = np.fft.fftshift(z)
    z = np.fft.fft2(z)
    z = np.fft.fftshift(z)
    z = z * z.conj()
    zmax = z.max()
    z = z/zmax

    return z

