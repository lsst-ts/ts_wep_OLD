import numpy as np
def getOffAxisCorr(confFile,fldr):

    #print confFile
    c=np.loadtxt(confFile)

    ruler = np.sqrt(c[:,0]**2+c[:,1]**2)
#    print ruler, fldr, (ruler >= fldr).argmax(), (ruler >= fldr).argmin()
    step=ruler[1]-ruler[0]

    p2=(ruler >= fldr)
#    print "FINE",p2, p2.shape
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
        p1=p2.argmax() - 1;
        p2 = p2.argmax()
        w1=(ruler[p2]-fldr)/step
        w2=(fldr-ruler[p1])/step

#    print c[p1,2:], c[p2,2:], w1,p1, w2, p2
    corr_coeff=np.dot(w1,c[p1,2:]) + np.dot(w2,c[p2,2:])

    return corr_coeff

