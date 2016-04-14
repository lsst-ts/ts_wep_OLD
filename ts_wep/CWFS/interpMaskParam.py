import numpy as np

def interpMaskParam(paramFile,fldr):
    
    c = np.loadtxt(paramFile)
    ruler = np.sqrt(2*c[:,0]**2)
    step=ruler[1]-ruler[0]

    p2=(ruler >= fldr)
    if (np.count_nonzero(p2) == 0):   #fldr is too large to be in range
        p2=c.shape()[0]
        p1=p2
        w1=1
        w2=0
    elif (p2[0]):   #fldr is too small to be in range
        p2 = 0 #this is going to be used as index
        p1=0 #this is going to be used as index
        w1=1
        w2=0
    else:
        p1=p2.argmax()-1
        p2 = p2.argmax()
        w1=(ruler[p2]-fldr)/step
        w2=(fldr-ruler[p1])/step

    param=np.dot(w1,c[p1,1:]) + np.dot(w2,c[p2,1:])
    ca = param[0]
    ra = param[1]
    cb = param[2]
    rb = param[3]

    return ca, ra, cb, rb

