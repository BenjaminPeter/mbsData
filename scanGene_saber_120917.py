#!/usr/bin/python

import sys
from mbsData import mbsData
import numpy as np

#1. argument is the mbs file you want to polarize
#2. argument is the outputfile in mbs format
#3. argument is the polariztion input in the following format:
# 1. col: chromosome
# 2. col: position
# 3. ool: human allele
# 4. col: chimp allele
# 5. col: maq allele
#4. set selected site

print len(sys.argv)
if len(sys.argv)<6:
    raise ValueError("""./scanGene_improved_10831.py [beagle output] [output file name] [file with polarizing info] [file with capture region info] [selected site] [mbs output file name]""")
windowSize=100000 #how large windows are
offset = 100000 #how much windows are apart
requirePolarizedData=True
#captReg=sys.argv[4]
#captReg="/home/peterb/Dropbox/mbsData/transfer/ABC_modifications/CaptReg.txt"
m = mbsData( )
m.readMsFile(sys.argv[1])
m.parTheta,m.parRho,m.nDataSets=0,0,0
#m.readRefSeq(sys.argv[3])
#m.polarizeData()
#m.setSequencedRegionTrack(captReg)

m.removeMonomorphicSites()
#m.removeSNPNotSequenced()
m.setSlidingWindowSets(windowSize=windowSize,offset=offset)

try:
	m.setSelectedSite(pos=int(sys.argv[5]))
except Error:
	print "error setting selected site"
	pass

m.writeToMbsFile(sys.argv[6])
#m.writeToMbsFile("%s.mbs"% (sys.argv[1]))

#toKeep=np.logical_and(m.segSites>lower,m.segSites<upper)

#m.data=m.data[:,toKeep]
#m.segSites=m.segSites[toKeep]
#m.nSegSites=sum(toKeep)

#m.createSNPdict()

td=m.getTajimasDfromSets()
print "1"
S=m.getSFromSets()
print "2"
pi=m.getPiFromSets()
print "3"
H=m.getHFromSets()
print "4"
EHH=m.getEHHforSetsAvg(windowSize)
print "5"
fwh=m.getFayWuHFromSets(requirePolarizedData=requirePolarizedData)
print "6"
#nsl=m.getNSLforSetsAvg(requirePolarizedData=requirePolarizedData)
nsl=np.zeros(len(m.sets))
print "7"
#IHS=m.getIHSforSetsAvg(requirePolarizedData=requirePolarizedData)
IHS=nsl
print "8"


output=np.empty((len(m.sets)+1,12),dtype="S20")
output[0:]=("startPos","endPos","midPos","coverage","td","S","pi","H","FWH","IHS","EHH","nsl")
output[1:,0]=m.windows[0]
output[1:,1]=m.windows[1]
output[1:,2]=m.windows[2]
output[1:,3]=m.coverage/float(m.windowSize)
output[1:,4]=td
output[1:,5]=S
output[1:,6]=pi
output[1:,7]=H
output[1:,8]=fwh
output[1:,9]=IHS
output[1:,10]=EHH
output[1:,11]=nsl
np.savetxt(sys.argv[2],output,fmt="%s",delimiter="\t")
