import numpy as np
def showProjection(lutxp, lutyp, sensorFactor, projSamples,raytrace):
    n1, n2 = lutxp.shape
    show_lutxyp = np.zeros((n1,n2))
    idx = (~np.isnan(lutxp)).nonzero()
#    idx = (~np.isnan(lutxp))
    for i,j in zip(idx[0],idx[1]):
        xR = np.round((lutxp[i,j] + sensorFactor) * \
                 (projSamples/sensorFactor)/2+0.5) #x=0.5 is center of pixel#1
        yR = np.round((lutyp[i,j] + sensorFactor) * (projSamples/sensorFactor)/2+0.5)
        
        if (xR>0 and xR<n2 and yR >0 and yR< n1):
            if raytrace:
                show_lutxyp[yR-1, xR-1] = show_lutxyp[yR-1, xR-1] +1
            else:
                if show_lutxyp[yR-1, xR-1] < 1:
                    show_lutxyp[yR-1, xR-1] = 1
    return show_lutxyp
