import os
import string

def readWavelength(wFile):
    """Read ..."""

    fid = open(os.path.join('wcs/conf/',wFile))

    wavelength = None
    for line in fid:
        if ((line.startswith('#') == False)  and len(line) > 0):
            if (line.startswith('PoissonSolver')):
                wavelength = string.atof(line.split()[1])
    fid.close()

    return wavelength
