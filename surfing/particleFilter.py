#!/usr/bin/env python
import numpy as np
import numpy.random as rng


def randomizePositions(nPositions,xlim=(0,100),ylim=(0,100)):
    return rng.uniform((xlim[0],ylim[0]),(xlim[1],ylim[1]),(nPositions,2))
    
def setJitteringConstants(Cr,Cq,k,R):
    Qbar=Cq/k/k
    Rbar=R+Cr/k/k
    return Qbar,Rbar

def h(x,y,D):
    deltaD=np.sqrt(y*y+(x+D/2)*(x+D/2))
    deltaD-=np.sqrt(y*y+(x-D/2)*(x-D/2))
    return deltaD

def hyperbolicFunctionGlobal(X,Y,Xi,Xj,Yi,Yj):
    alpha=np.arctan((Yi-Yj)/(Xi-Xj))
    X0=(Xi+Xj)/2
    Y0=(Yi-Yj)/2
    x=(X-X0)*np.cos(alpha)+(Y-Y0)*np.sin(alpha)
    y=-1*(X-X0)*np.sin(alpha)+(Y-Y0)*np.cos(alpha)
    D=np.sqrt((Yi-Yj)**2+(Xi-Xj)**2)
    print alpha,X0,Y0,x,y,D
    return x,y,D,h(x,y,D)
                    

def makeHyperbolaFunction(s1,s2,measurement,variance):
    def lsqerror(x,y,v):
        error=0.0
        for i in range(len(measurement)):
            error+=variance[i](v*measurement[i]-hyperbolicFunctionGlobal(x[i],y[i],
                                        s1[i,0],s2[i,0],s1[i,1],s2[i,1]))**2
        return error





def computeParticleWeights(particles):
    logweights=error

def normalize(particleWeights):
    particleWeights/=sum(particleWeights)

def resample(positions,particleWeights):
    which=rng.multinomial(len(positions),particleWeights)
    return np.repeat(positions,which)

def spreadOut(positions,Qbar):
    noise=rng.normal(positions,Qbar)


def doIterations(kMax,positions,Cr,Cq):
    for k in range(kMax):
        Qbar,Rbar=setJitteringConstants(Cr,Cq,k,R)
        particleWeights=computeParticleWeights(particles, weightFunction) #step a.1
        particleWeights=normalize(particleWeights) #step a.2
        estimate=computeEstimate(particleWeights,positions) #step b
        positions=resample(positions,particleWeights) #step c
        positions=spreadOut(positions,Qbar) #step d
    return estimate

#------------------------
nPositions=100
kMax=2
positions=randomizePositions(nPositions) #step 1
Cr,Cq = 1, 1 #step 2 constants=setJitteringConstants(Cr,Cq) #step 2

result=doIterations(kMax,positions,Cr,Cq) #step 3

