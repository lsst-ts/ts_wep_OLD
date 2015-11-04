import numpy as np
def upResolution(myI,oversample,lm,ln):
    
    #lm and ln are dimensions after upResolution
    sm=lm/oversample
    sn=ln/oversample

    newI=np.zeros(lm,ln);
    for i in range(sm):
        for j in range(sn):
            for k in range(oversample):
                for l in range(oversample):
                    newI[(i-1)*oversample+k,(j-1)*oversample+l]=myI[i,j]/oversample/oversample;
    return newI
