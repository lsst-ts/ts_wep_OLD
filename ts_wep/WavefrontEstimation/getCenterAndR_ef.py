import numpy as np
import string

def getCenterAndR_ef(oriArray,readRand=1):
    #centering finding code based northcott_ef_bundle/ef/ef/efimageFunc.cc
    #this is the modified version of getCenterAndR_ef.m 6/25/14

    stepsize = 20
    nwalk = 400
    slide=220
    
    histogram_len = 256

    array=oriArray.copy()
    array[array<0]=0; #deblending can make over-subtraction in some area, a lower tail could form; temperary solution.
    m, n = array.shape

    pmin = array.min()
    pmax = array.max()
    if (pmin == pmax):
        print('image has min=max=%f',pmin)
    array1d = np.reshape(array,m*n,1)

    phist, cen = np.histogram(array1d,bins=histogram_len)
    startidx=range(25,175+25,25)
    # startidx=fliplr(startidx)
    foundvalley=0
    
    #validating code against Matlab:
    #to avoid differences in random number generator between NumPy and Matlab, read in these random numbers generated from Matlab
    if readRand:
        iRand=0
        myRand = np.loadtxt('../conf/testRand.txt')
        myRand=np.tile(np.reshape(myRand,(1000,1)),(10,1))

    for istartPoint in range(len(startidx)):
        minind = startidx[istartPoint]
        if ((minind <=0) or (max(phist[minind-1:])==0)):
            continue
        minval = phist[minind-1]

        #try nwalk=2000 times and see if it rolls to > slide. 
        #if it does roll off, let's change starting point (minind) to 25 less
        #(minind starts at 175, then goes to 150, then 125
        for i in range(nwalk+1):
            if (minind <= slide):
                while (minval==0): 
                    minind=minind+1
                    minval=phist[minind-1]
                if readRand:
                     ind = np.round(-stepsize+2*stepsize*myRand[iRand,0]);
                     iRand += 1
                     thermal = 1+0.5*myRand[iRand,0]*np.exp(-(1.0*i/(nwalk*0.3)))
                     iRand += 1
                else:
                     ind = np.round(stepsize*(2*np.random.rand()-1))
                     thermal = 1 #+0.05*np.random.rand() #*np.exp(-(1.0*i/(nwalk*0.3)))
                if ((minind+ind<1) or (minind+ind > (histogram_len))):
                    continue
                if (phist[minind+ind-1] < (minval*thermal)):
                    minval = phist[minind+ind-1]
                    minind = minind + ind
            else:
                break
        if (minind <= slide):
            foundvalley+=1
        if foundvalley>=1:
            break

    # Never understood why we do this, but had it before 5/27/14, since EF C++
    # code had this. Now this appears to cause minind> histgram_length, we
    # comment it out and see how the code performs.
    # minind = avind/(steps-steps/2);

    # fprintf('minind=%d\n',minind);
    if (not foundvalley): #because of noise, there is only peak
        minind = histogram_len/2
    pval = pmin+(pmax-pmin)/histogram_len*float(minind-0.5)
    tmp = array.copy()
    tmp[array>max(0,pval-1e-8)] = 1
    tmp[array<pval] = 0
    realR = np.sqrt(np.sum(tmp)/3.1415926) #because tmp is a binary mask with only 1 and 0

    jarray, iarray = np.mgrid[1:n+1,1:m+1]
    realcx=np.sum(iarray*tmp)/np.sum(tmp)
    realcy=np.sum(jarray*tmp)/np.sum(tmp)

    #print realcx, realcy, realR
    return realcx, realcy, realR
