#!/usr/bin/env python

from mbsData import mbsData
from mbsDataMultipop import mbsDataMP
import numpy as np

class mbsDataMPSSS(mbsDataMP):
    """this should allow usage of the standard mbsData function using the output from the SSSimulator. Currently, data comes in the following form:
        1. 2D SFS, ending *.sfs
            this files have a 2dSFS, the first line is the total tree length, the second through nth lines are the actual SFS
        2. Frequency tables, ending *.ft
            each row of a FT is a tree segment, giving the lenght of the tree segment followed by the (absolute) allele frequency in each population. Sample Sizes have to be given externally...
        3. output of statistics, ending *.stat
            gives popid1, popid2 FST sdFST, #psi sdpsi, deltaH, sddeltahfor all pairwise comparisons

10.9.2012: currently solely working on FT
    """
#--------------------------------------------------
    def readFT(self,file="fig3.e/e100_67.ft",ssfile="sampSizes.ss"):
        """reads file format where the first column is the branch length, following columns are the absolute frequency in all populations. Currently sample sizes in all pops are assumed to be 100
        in contrast to the mbsDataMultipop, data is stored as allele frequency of each branch in self.data, and the length of the branches .in self.dataFreq"""
        self.data=np.loadtxt(file)
        self.dataFreq=self.data[:,0]
        self.data=self.data[:,1:]

        #nSegsites here is just the number of entries
        self.nSegsites=self.data.shape[0] 
        self.tTot=sum(self.dataFreq)
        self.nPops=self.data.shape[1]#-1 due to first col being branch length

        self.ss=np.zeros(self.nPops-1,dtype="i4")+100
        self.setActivePops([0,1])
#--------------------------------------------------
    def setPops(self, pops):
        """pops should be a list of lists of the IDs for each individual in the population, overrides function in mbsDataMultipop.py"""
        self.nPops=len(pops)
        self.pops=[np.array(p) for p in pops]
        self.popVector=np.zeros(self.nHap,dtype="i4")
        for i,p in enumerate(pops):
            for pp in p:
                self.popVector[pp]=i+1
        self.n1=len(self.pops[0])
        self.n2=len(self.pops[1])
#--------------------------------------------------
    def setActivePops(self,pops,resample=False):
        """overrides the corresponding function in mbsDataMultipop.py"""
        i,j=tuple(pops)
        self.activePops=pops
        self.nHap=self.ss[i]+self.ss[j]
        self.setPops((range(self.ss[i]),range(self.ss[i],self.ss[i]+self.ss[j])))
        if hasattr(self,"data"):
            self.createMultiDimSFS(resample=resample)
        if hasattr(self,"sets"):
            tmpSets=self.sets
            self.sets=self.blocks
            self.createMDSFSForSets()
            for i in range(len(self.setMDSFS)):
                self.setMDSFS[i]=self.mdsfs-self.setMDSFS[i]
            self.sets=tmpSets
        if resample:
            self.n1 = min(self.n1,self.n2)
            self.n2 = self.n1
#--------------------------------------------------
    def createMultiDimSFS(self,resample=False):
        """creates 2dsfs for 2 populations, if resample is set to true, the larger population will be downsampled so that sample sizes are equal, override"""
        self.mdfreq=np.transpose(self.data[:,(self.activePops)])
        if resample and self.n1 != self.n2:
            if self.n1>self.n2:
                self.mdfreq[0]=[np.random.binomial(self.n2,float(f)/self.n1) \
                                for f in self.mdfreq[0]]
            else:
                self.mdfreq[1]=[np.random.binomial(self.n1,float(f)/self.n2) \
                                for f in self.mdfreq[1]]

        mdsfs=np.zeros((self.ss[self.activePops[0]]+1,\
                        self.ss[self.activePops[1]]+1),\
                        dtype="f4")

        shared=0.
        private1=0.
        private2=0.
        for i in range(self.nSegsites):
            print tuple(self.mdfreq[:,i])
            mdsfs[tuple(self.mdfreq[:,i])]+=self.dataFreq[i]
            if len(mdsfs.shape)==2:
                if self.mdfreq[0,i]>0 and self.mdfreq[1,i]>0:
                    shared+=self.dataFreq[i]
                elif self.mdfreq[0,i]>0:
                    private1+=self.dataFreq[i]
                elif self.mdfreq[1,i]>0:
                    private2+=self.dataFreq[i]
                else:
                    pass
                    #print self.mdfreq[:,i],
                    #print "ERROR: neither private nor shared "

        self.totNAlleles=sum(sum(mdsfs))
        #self.shared=shared/self.totNAlleles
        #self.private=(private1/self.totNAlleles,private2/self.totNAlleles)
        #self.mdsfs=mdsfs/self.totNAlleles
        return mdsfs
#--------------------------------------------------
    def createMultiDimSFSFolded(self,resample=False):
        """creates 2dsfs for 2 populations, if resample is set to true, the larger population will be downsampled so that sample sizes are equal, override"""
        self.mdfreq=np.transpose(self.data[:,(self.activePops)])
        if resample and self.n1 != self.n2:
            if self.n1>self.n2:
                self.mdfreq[0]=[np.random.binomial(self.n2,float(f)/self.n1) \
                                for f in self.mdfreq[0]]
            else:
                self.mdfreq[1]=[np.random.binomial(self.n1,float(f)/self.n2) \
                                for f in self.mdfreq[1]]

        mdsfs=np.zeros((self.ss[self.activePops[0]]+1,\
                        self.ss[self.activePops[1]]+1),\
                        dtype="f4")

        shared=0.
        private1=0.
        private2=0.
        for i in range(self.nSegsites):
            mdsfs[tuple(self.mdfreq[:,i])]+=self.dataFreq[i]
            if len(mdsfs.shape)==2:
                if self.mdfreq[0,i]>0 and self.mdfreq[1,i]>0:
                    shared+=self.dataFreq[i]
                elif self.mdfreq[0,i]>0:
                    private1+=self.dataFreq[i]
                elif self.mdfreq[1,i]>0:
                    private2+=self.dataFreq[i]
                else:
                    pass
                    #print self.mdfreq[:,i],
                    #print "ERROR: neither private nor shared "

        self.totNAlleles=sum(sum(mdsfs))
        #self.shared=shared/self.totNAlleles
        #self.private=(private1/self.totNAlleles,private2/self.totNAlleles)
        #self.mdsfs=mdsfs/self.totNAlleles
        return mdsfs
#--------------------------------------------------


    def createMDSFSForSets(self,folded=False):
        self.setMDSFS=[]
        for s in self.sets:
            sfs=np.zeros((self.n1+1,self.n2+1),dtype="int")
            for i in range(self.nSegSites):
                if i in s:
                    sfs[tuple(self.mdfreq[:,i])]+=self.dataFreq[i]
            self.setMDSFS.append(sfs)
        if self.setsIncludeGlobal:
            self.setSFS.insert(0,self.mdsfs)
        self.setMDSFS=np.array(self.setMDSFS)

    def createSFSForSets(self):
        self.setSFS=[]
        for s in self.sets:
            sfs=[np.zeros(self.n1+1),
                 np.zeros(self.n2+1)
                 np.zeros(self.nHap)]
            for i in np.arange(self.nSegSites):
                if i in s:
                    sfs[0][self.mdfreq[0,i]]+=self.dataFreq[i]
                    sfs[1][self.mdfreq[0,i]]+=self.dataFreq[i]
                    sfs[2][sum(self.mdfreq[:,i]])+=self.dataFreq[i]
            self.setSFS.append(sfs)
        if self.setsIncludeGlobal:
            sfs=[np.zeros(self.n1+1),
                 np.zeros(self.n2+1)
                 np.zeros(self.nHap)]
            for i in np.arange(self.nSegSites):
                sfs[0][self.mdfreq[0,i]]+=self.dataFreq[i]
                sfs[1][self.mdfreq[0,i]]+=self.dataFreq[i]
                sfs[2][sum(self.mdfreq[:,i]])+=self.dataFreq[i]
            self.setSFS.append(sfs)
