import numpy

def downResolution(myI,oversample,sm,sn):
#sm and sn are dimensions after downResolution

    newI= numpy.Zeros(sm,sn)
    for i in range(1,sm)
        for j in range(1,sn):
            for k in range(1,oversample):
                for l in range(1,oversample):
                    newI[i,j]=newI[i,j]+myI[(i-1)*oversample+k,(j-1)*oversample+l]
    
    return newI
