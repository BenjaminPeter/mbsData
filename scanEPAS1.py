#!/usr/bin/python
from mbsData import mbsData
import numpy as np
import sys

if len(sys.argv)<4:
    raise ValueError("call: ./this dataset(mbs format), output filter(0=no,1=yes)")
tib = mbsData(sys.argv[1])
#tib.readPhasedFile("phased/Epas1_v3.tibet.like.phased",filter=False)
#tib.nHap=len(tib.tib)
#tib.data=tib.data[tib.tib,:]

print len(tib.segSites)
print tib.nSegsites
print tib.data.shape
if sys.argv[3]=="1":
    tib.removeSingletons()
    tib.removeDoubletons()
"""tib.sfs,tib.freq=tib.createSFS()
toRemove=np.where(tib.freq>2)[0]
tib.data=tib.data[:,toRemove]
tib.segSites=tib.segSites[toRemove]
tib.nSegsites=len(tib.segSites)

tib.createSNPdict()
tib.selId=tib.SNPdict[tib.selPos]
tib.selSiteData=tib.selSiteData[tib.tib]"""

windowSize=10000
td=np.empty(tib.nSegsites)
S=np.empty(tib.nSegsites)
Pi=np.empty(tib.nSegsites)
H=np.empty(tib.nSegsites)
IHS=np.empty(tib.nSegsites)
rg=np.empty(tib.nSegsites)
EHH=np.empty(tib.nSegsites)
uHap=np.empty(tib.nSegsites)
fwh=np.empty(tib.nSegsites)
daf=np.empty(tib.nSegsites)
#k                       =       tib.getNSL()
#sl=np.log(tib.nsl[1]/tib.nsl[2])

for snp in np.arange(tib.nSegsites):
    pos                     =       tib.segSites[snp]
    interval                =       (pos-windowSize/2,pos+windowSize/2)
    ind                     =       tib.getIndividuals()
    td[snp],S[snp],Pi[snp]  =       tib.getTajimasD(individuals=ind,pos=interval)
    H[snp]                  =       tib.getH(pos=interval)
    #IHS[snp],_,_,_,(dmin,dmax)    =       tib.getIHS(id=snp)
    IHS[snp]                =       0
    EHH[snp]                =       tib.getEHH(windowSize/2,id=snp)
    hap                     =       tib._getUniqueHaplotypes(pos=interval).values()
    uHap[snp]               =       len(hap)
    rg[snp]              =       np.max(hap)
    fwh[snp]             =          tib.getFayWuH(pos=interval)
    daf[snp]                =       len(tib.getIndividualsWithDerivedAllel(id=snp))
    print snp,pos,interval,daf[snp]
    
#masonStatistic(sys.argv[1],"temp")
#sl=np.loadtxt("temp")

output=np.empty((tib.nSegsites+1,11),dtype="S20")
output[0:]=("pos","td","S","pi","H","FWH","IHS","EHH","uHap","maxHap","daf")
output[1:,0]=tib.segSites
output[1:,1]=td
output[1:,2]=S
output[1:,3]=Pi
output[1:,4]=H
output[1:,5]=fwh
output[1:,6]=IHS
output[1:,7]=EHH
output[1:,8]=uHap
output[1:,9]=rg
output[1:,10]=daf

np.savetxt(sys.argv[2],output,fmt="%s",delimiter="\t")


