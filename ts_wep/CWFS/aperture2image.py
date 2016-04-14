import numpy as np
from ZernikeAnnularGrad import ZernikeAnnularGrad
from ZernikeGrad import ZernikeGrad
from ZernikeJacobian import ZernikeJacobian
from ZernikeAnnularJacobian import ZernikeAnnularJacobian
from createPupilGrid import createPupilGrid

def aperture2image(m,zcCol,lutx,luty,projSamples,atype,fldx,fldy,model):

    R = m.apertureDiameter/2.0
    if (atype == 'intra'):
        myC = - m.focalLength*(m.focalLength-m.offset)/m.offset/R**2
    elif (atype == 'extra'):
        myC =  m.focalLength*(m.focalLength/m.offset+1)/R**2

    module = __import__ ('CWFS.poly%d_2D'%m.offAxisPolyOrder)
    polyFunc = getattr(module, 'poly%d_2D'%m.offAxisPolyOrder)

    module = __import__ ('CWFS.poly%dGrad'%m.offAxisPolyOrder)
    polyGradFunc = getattr(module, 'poly%dGrad'%m.offAxisPolyOrder)

    lutr=np.sqrt(lutx**2+luty**2)
    #1 pixel larger than projected pupil. No need to be EF-like, anything
    #outside of this will be masked off by the computational mask
    onepixel=1/(projSamples/2/m.sensorFactor)
    #print "LUTR",lutr.shape,((lutr >1 + onepixel) | (lutr<m.obscuration-onepixel)).shape
    idxout=((lutr >1 + onepixel) | (lutr<m.obscuration-onepixel))
    lutx[idxout]=np.nan
    luty[idxout]=np.nan
    idxbound= ((lutr<=1+onepixel) & (lutr>1))#outer boundary (1 pixel wide boundary)
    lutx[idxbound]=lutx[idxbound]/lutr[idxbound]
    luty[idxbound]=luty[idxbound]/lutr[idxbound]
    idxinbd=((lutr<m.obscuration) & (lutr>m.obscuration-onepixel)) #inner boundary
    lutx[idxinbd] = lutx[idxinbd]/lutr[idxinbd]*m.obscuration
    luty[idxinbd] = luty[idxinbd]/lutr[idxinbd]*m.obscuration

    if (model == 'offAxis'):
        if (atype == 'intra'):
            ca=m.I1maskCa
            cb=m.I1maskCb
            ra=m.I1maskRa
            rb=m.I1maskRb
        elif (atype == 'extra'):
            ca=m.I2maskCa
            cb=m.I2maskCb
            ra=m.I2maskRa
            rb=m.I2maskRb

        lutx, luty = createPupilGrid(lutx,luty,onepixel,ca,cb,ra,rb,fldx,fldy)

    if (model == 'paraxial'):
        lutxp = lutx
        lutyp = luty
    elif (model == 'onAxis'):
        myA=np.sqrt((m.focalLength**2-R**2)/(m.focalLength**2-lutr**2*R**2))
        lutxp = m.maskScalingFactor*myA*lutx
        lutyp = m.maskScalingFactor*myA*luty
    elif (model=='offAxis'):
        #nothing is hard-coded here: (1e-3) is because the
        # offAxis correction param files are based on offset=1.0mm
        if (atype == 'intra'):
            cx=(m.I1offAxis_coeff[0,:]-m.I1offAxis_coeff[2,:])*(1e-3+m.offset)/2e-3+m.I1offAxis_coeff[2,:]
            cy=(m.I1offAxis_coeff[1,:]-m.I1offAxis_coeff[3,:])*(1e-3+m.offset)/2e-3+m.I1offAxis_coeff[3,:]
        elif (atype == 'extra'):
            cx=(m.I2offAxis_coeff[0,:]-m.I2offAxis_coeff[2,:])*(1e-3-m.offset)/2e-3+m.I2offAxis_coeff[2,:]
            cy=(m.I2offAxis_coeff[1,:]-m.I2offAxis_coeff[3,:])*(1e-3-m.offset)/2e-3+m.I2offAxis_coeff[3,:]
            cx=-cx #this will be inverted back by typesign later on.
            cy=-cy #we do the inversion here to make the (x,y)->(x',y') equations has the same form as the paraxial case.

        fldr = np.sqrt(fldx**2+fldy**2)
        costheta = (fldx+fldy)/fldr/1.4142
        if (costheta>1):
            costheta=1
        elif (costheta<-1):
            costheta = -1

        sintheta = np.sqrt(1-costheta**2)
        if (fldy<fldx):
            sintheta = -sintheta

        lutx0 = lutx*costheta + luty*sintheta #first rotate back to reference orientation
        luty0 = -lutx*sintheta + luty*costheta
        lutxp0 = polyFunc(cx,lutx0,y=luty0) #use mapping at reference orientation
        lutyp0 = polyFunc(cy,lutx0,y=luty0)
        lutxp = lutxp0*costheta - lutyp0*sintheta #rotate back
        lutyp = lutxp0*sintheta + lutyp0*costheta
        reduced_coordi_factor = 1/(m.sensorSamples/2*m.pixelSize/m.sensorFactor)/1000; #Zemax data are in mm, therefore 1000
        lutxp=lutxp*reduced_coordi_factor #reduced coordinates, so that this can be added with the dW/dz
        lutyp=lutyp*reduced_coordi_factor
    else:
        printf('wrong model number in compensate\n')
        return


    if (zcCol.ndim ==1):
        if (m.zobsR>0):
            lutxp =lutxp + myC* ZernikeAnnularGrad(zcCol,lutx, luty,m.zobsR,'dx')
            lutyp =lutyp + myC* ZernikeAnnularGrad(zcCol,lutx, luty,m.zobsR,'dy')
        else:
            lutxp = lutxp + myC * ZernikeGrad(zcCol,lutx, luty,'dx')
            lutyp = lutyp + myC * ZernikeGrad(zcCol,lutx, luty,'dy')
    else:
        FX, FY=gradient(zcCol,m.sensorFactor/(m.sensorSamples/2))
        lutxp =lutxp + myC* FX
        lutyp =lutyp + myC* FY


    if  (atype == 'extra'):
        lutxp = - lutxp
        lutyp = - lutyp


    #%%%%Below for calculation of the Jacobian

    if (zcCol.ndim ==1):
        if (model == 'paraxial'):
            if (m.zobsR>0):
                J= (1+ myC * ZernikeAnnularJacobian(zcCol, lutx, luty,m.zobsR, '1st') +
                    myC**2 * ZernikeAnnularJacobian(zcCol, lutx, luty,m.zobsR, '2nd'))
            else:
                J= (1+ myC * ZernikeJacobian(zcCol, lutx, luty, '1st') +
                    myC**2 * ZernikeJacobian(zcCol, lutx, luty, '2nd'))

        elif (model == 'onAxis'):
            xpox = m.maskScalingFactor*myA*(1+lutx**2*R**2./(m.focalLength**2-R**2*lutr**2)) + \
                   myC * ZernikeAnnularGrad(zcCol,lutx, luty,m.zobsR,'dx2')
            ypoy = m.maskScalingFactor*myA*(1+luty**2*R**2./(m.focalLength**2-R**2*lutr**2)) \
                   +myC* ZernikeAnnularGrad(zcCol,lutx, luty,m.zobsR,'dy2')
            xpoy = m.maskScalingFactor*myA*lutx*luty*R**2/(m.focalLength**2-R**2*lutr**2) \
                   +myC*ZernikeAnnularGrad(zcCol,lutx, luty,m.zobsR,'dxy')
            ypox = xpoy
            J = ( xpox*ypoy - xpoy*ypox)
        elif (model == 'offAxis'):
            xp0ox=polyGradFunc(cx,lutx0,luty0,'dx')*costheta-polyGradFunc(cx,lutx0,luty0,'dy')*sintheta
            yp0ox=polyGradFunc(cy,lutx0,luty0,'dx')*costheta-polyGradFunc(cy,lutx0,luty0,'dy')*sintheta
            xp0oy=polyGradFunc(cx,lutx0,luty0,'dx')*sintheta+polyGradFunc(cx,lutx0,luty0,'dy')*costheta
            yp0oy=polyGradFunc(cy,lutx0,luty0,'dx')*sintheta+polyGradFunc(cy,lutx0,luty0,'dy')*costheta
            xpox = (xp0ox*costheta-yp0ox*sintheta) *reduced_coordi_factor \
                   +myC* ZernikeAnnularGrad(zcCol,lutx, luty,m.zobsR,'dx2')
            ypoy = (xp0oy*sintheta+yp0oy*costheta)*reduced_coordi_factor \
                   +myC* ZernikeAnnularGrad(zcCol,lutx, luty,m.zobsR,'dy2')
            temp=myC*ZernikeAnnularGrad(zcCol,lutx, luty,m.zobsR,'dxy')
            xpoy = (xp0oy*costheta-yp0oy*sintheta)*reduced_coordi_factor+temp #if temp==0,xpoy doesn't need to be symmetric about x=y
            ypox = (xp0ox*sintheta+yp0ox*costheta)*reduced_coordi_factor+temp #xpoy-flipud(rot90(ypox))==0 is true
            J= ( xpox*ypoy - xpoy*ypox)

    else:

        FXX, FXY =gradient(FX,m.sensorFactor/(m.sensorSamples/2))
        tmp, FYY =gradient(FY,m.sensorFactor/(m.sensorSamples/2))
        if (model == 'paraxial'):
            xpox = 1+myC* FXX
            ypoy = 1+myC* FYY
            xpoy = 1+myC* FXY
            ypox = xpoy
        elif (model == 'onAxis'):
            xpox = m.maskScalingFactor*myA*(1+lutx**2*R**2./(m.focalLength**2-R**2*lutr**2)) +myC* FXX
            ypoy = m.maskScalingFactor*myA*(1+luty**2*R**2/(m.focalLength**2-R**2*lutr**2)) +myC* FYY
            xpoy = m.maskScalingFactor*myA*lutx*luty*R**2./(m.focalLength**2-R**2*lutr**2)   +myC* FXY
            ypox = xpoy
        elif (model == 'offAxis'):
            xpox = polyGradFunc(cx,lutx,luty,'dx')*reduced_coordi_factor +myC* FXX
            ypoy = polyGradFunc(cy,lutx,luty,'dy')*reduced_coordi_factor  +myC* FYY
            xpoy = polyGradFunc(cx,lutx,luty,'dy')*reduced_coordi_factor  +myC* FXY
            ypox = polyGradFunc(cy,lutx,luty,'dx')*reduced_coordi_factor  +myC* FXY

        J= ( xpox*ypoy - xpoy*ypox)

    return lutxp, lutyp, J
