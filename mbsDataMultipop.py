#!/usr/bin/python

from mbsData import mbsData
import numpy as np
import os
import numpy.random as rng

class mbsDataMP(mbsData):
    def __init__(self):
        self.nPops=1
        self.pops       =   None
        self.popVector  =   None

    def readSplatche(self,file="sGeneticsOutput/_spc_samp_1.arp",
                     samfile="sc_samp.sam"):
        self.readArpFile(file)
        f=open(samfile)
        nSamplesExp=int(f.readline().split()[0])
        self.pos=np.zeros((nSamplesExp,2),dtype="i")
        expectedSamples=np.zeros(nSamplesExp,dtype="S20")
        for i in range(nSamplesExp):
            line=f.readline()
            self.pos[i]=line.split()[3:6]
            expectedSamples[i]=line.split()[0]
        f.close()

        foundSamples=np.array([np.where(i==expectedSamples)[0] for i in
                              self.sampleId]).flatten()
        #remove samples that are not present
        self.pos=self.pos[foundSamples]
        self.nPops=len(self.ss)

        data=self.data
        self.data=np.empty((self.nPops,self.nSegsites))
        cummulativePos=0
        for i,popSize in enumerate(self.ss):
            self.data[i]    =   sum(data[range(cummulativePos,
                                               cummulativePos+popSize)])
            cummulativePos+=popSize
        self.data=np.transpose(self.data)


    def getPairwiseStat(self,stat):
        statMat=np.empty((self.nPops,self.nPops))
        for i in range(self.nPops):
            for j in range(self.nPops):
                self.setActivePops((i,j))
                statMat[i,j]=stat()
        return statMat

        
    def readHGDP(self,f1="/data/surfing/applications/hgdp/metric/daf_cnt.txt",
                f2="/data/surfing/applications/hgdp/metric/popSS.txt"):
        self.data=np.loadtxt(f1,dtype="i4",skiprows=000000)
        self.popData=np.loadtxt(f2,dtype="S")
        self.ss=2*np.array([int(i) for i in self.popData[:,1]])
        self.nSegsites=self.data.shape[0]
        self.setActivePops([0,1])
        
    def setPops(self, pops):
        """pops should be a list of lists of the IDs for each individual in the population"""
        self.nPops=len(pops)
        self.pops=[np.array(p) for p in pops]
        self.popVector=np.zeros(self.nHap,dtype="i4")
        for i,p in enumerate(pops):
            for pp in p:
                self.popVector[pp]=i+1
        self.n1=len(self.pops[0])
        self.n2=len(self.pops[1])
        #print self.popVector
    
    def setActivePops(self,pops,resample=False):
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

    def createMultiDimSFS(self,resample=False):
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
                        dtype="I")

        shared=0
        private1=0
        private2=0
        for i in range(self.nSegsites):
            mdsfs[tuple(self.mdfreq[:,i])]+=1
            if len(mdsfs.shape)==2:
                if self.mdfreq[0,i]>0 and self.mdfreq[1,i]>0:
                    shared+=1
                elif self.mdfreq[0,i]>0:
                    private1+=1
                elif self.mdfreq[1,i]>0:
                    private2+=1
                else:
                    pass
                    #print self.mdfreq[:,i],
                    #print "ERROR: neither private nor shared "

        self.shared=shared
        self.private=(private1,private2)
        self.mdsfs=mdsfs
        return mdsfs
    def createMultiDimSFSFolded(self,folded=False):
        mdfreq=np.empty((self.nPops,self.nSegsites),dtype="I")
        for i in range(self.nSegsites):
            snp=self.getDataFromId(id=i)
            if folded:
                if sum(snp)>sum(1-snp):
#                    print "inverted:[%d -> %d]"%(sum(snp),sum(1-snp))
                    snp=1-snp
            for j in range(self.nPops):
                #print i,j,snp
                mdfreq[j,i]=sum(snp[self.pops[j]])

        dimensions=[len(p)+1 for p in self.pops]
        #print dimensions
        mdsfs=np.zeros(dimensions,dtype="I")

        shared=0
        private1=0
        private2=0
        for i in range(self.nSegsites):
            mdsfs[tuple(mdfreq[:,i])]+=1
            if len(mdsfs.shape)==2:
                if mdfreq[0,i]>0 and mdfreq[1,i]>0:
                    shared+=1
                elif mdfreq[0,i]>0:
                    private1+=1
                elif mdfreq[1,i]>0:
                    private2+=1
                else:
                    print "ERROR: neither private nor shared "

        self.shared=shared
        self.private=(private1,private2)
        self.mdfreq,self.mdsfs=mdfreq,mdsfs
        return mdfreq,mdsfs

    def readHGDPBlocks(self,filter=None,file="/data/surfing/applications/hgdp/metric/snp2.txt"):
        sets=dict()
        f=open(file,"r")
        i=0
        while True:
            line=f.readline()
            if line=="":
                break
            chr,pos,id,s=line.split()
            if not int(s) in sets:
                sets[int(s)]=[]
            if filter==None:
                sets[int(s)].append(i)
            elif i==filter[retainedPos]:
                sets[int(s)].append(retainedPos)
                retainedPos+=1
            else:
                print "filtered ",i,line
            i+=1
        f.close()
        self.blocks=[np.array(s) for s in sets.values()]
        self.sets=[np.array(s) for s in sets.values()]
        self.setsIncludeGlobal=False
        self.createMultiDimSFS()
        self.createMDSFSForSets()
        for i in range(len(self.setMDSFS)):
            self.setMDSFS[i]=self.mdsfs-self.setMDSFS[i]
            
        for i in range(len(self.sets)):
            self.sets[i]=np.setdiff1d(np.arange(self.nSegsites),self.sets[i])

    def readNeurospora(self):
        m.data=np.loadtxt("bla")
        m.ss=[20,19]
        m.nSegsites=m.data.shape[0]
        m.setActivePops((0,1))
    def _getSurfStatData(self,sfs,l,u,b,shared,verbose):
        k=0.
        for i in range(l,u[0]):
            for j in range(l,u[1]):
                k+=sfs[i,j]*(b[0]*i-b[1]*j)
                if verbose: print k,i,j,sfs[i,j]
        return float(k)/shared

    def _getSurfStat(self,l=1,u=0,b=(1.,1.),verbose=True):
        if not hasattr(self,"mdsfs"):
            self.createMultiDimSFS()
        if self.shared==0:
            return 0

        return self._getSurfStatData(self.mdsfs,l,(self.n1+u,self.n2+u),b,self.shared,verbose)

    def getPhi1(self,verbose=False):
        return self._getSurfStat(1,0,(1.,1.),verbose)

    def getPhi2(self,verbose=False):
        n1=self.n1; n2=self.n2
        return self._getSurfStat(1,1,((1.+n2)/(n1*n2),(n1+1.)/(n1*n2)),verbose)

    def getPhi3(self,verbose=False):
        n1=self.n1; n2=self.n2
        return self._getSurfStat(1,1,(1./n1,1./n2),verbose)

    def getPhi4(self,verbose=False):
        n1=self.n1; n2=self.n2
        return self._getSurfStat(1,0,((1.+n2)/(n1*n2),(n1+1.)/(n1*n2)),verbose)

    def _getSurfStatForSets(self,l=1,u=0,b=(1.,1.),verbose=True):
        phi=[]
        for sfs in self.setMDSFS:
            shared=sum(sum(sfs[1:,1:]))
            phi.append(self._getSurfStatData(sfs,l,(self.n1+u,self.n2+u),b,shared,verbose))
        return phi

    def getPhi1ForSets(self,verbose=False):
        return self._getSurfStatForSets(l=1,u=0,b=(1.,1.),verbose=verbose)

    def getPhi2ForSets(self,verbose=False):
        n1=self.n1; n2=self.n2
        return self._getSurfStatForSets(1,1,((1.+n2)/(n1*n2),(n1+1.)/(n1*n2)),verbose=verbose)

    def getPhi3ForSets(self,verbose=False):
        n1=self.n1; n2=self.n2
        return self._getSurfStatForSets(0,0,(1./n1,1./n2),verbose=verbose)

    def getPhi4ForSets(self,verbose=False):
        n1=self.n1; n2=self.n2
        return self._getSurfStatForSets(1,0,((1.+n2)/(n1*n2),(n1+1.)/(n1*n2)),verbose=verbose)

    def getDeltaH(self):
        h1=self.getH(individuals=self.pops[0])
        h2=self.getH(individuals=self.pops[1])        
        return h2-h1
    def create1dsfsFrom2dsfs(self,sfs,n1,n2):
        sfs1=np.zeros(n1+1,dtype="i4")
        sfs2=np.zeros(n2+1,dtype="i4")
        S1,S2=0,0
        for i in range(n1+1):
            for j in range(n2+1):
                sfs1[i]+=sfs[i,j]
                sfs2[j]+=sfs[i,j]
                if i>0 and i < n1:
                    S1+=sfs[i,j]
                if j>0 and j < n2:
                    S2+=sfs[i,j]
        return sfs1,sfs2, S1,S2   
        
    def _getPi(self,sfs,n):           
        denominator = n*(n-1.0)

        thetaPi=0.0
        for i,freq in enumerate(sfs):
            if i==0 or i==n:
                continue
            thetaPi+= 2.0*i*(n-i)*freq
        return thetaPi/denominator
             
    def getDeltaHForSets(self):
    
        deltaH=[]
        for sfs in self.setMDSFS:
            sfs1,sfs2,S1,S2=self.create1dsfsFrom2dsfs(sfs,self.n1,self.n2)
            h1=self._getPi(sfs1,self.n1)/S1
            h2=self._getPi(sfs2,self.n2)/S2
#            h1=self.getH(id=s,individuals=self.pops[0])
#            h2=self.getH(id=s,individuals=self.pops[1])
            deltaH.append(h2-h1)
        return np.array(deltaH)

    def getFST(self):
        fst=0.0
        for i in range(self.n1+1):
            for j in range(self.n2+1):
                fst += self.mdsfs[i,j]*self.getFST1Locus(float(i)/self.n1,float(j)/self.n2,self.n1/2.,self.n2/2.)
        return fst/self.nSegsites

    def getFST1Locus(self,p1,p2,n1,n2):
        alpha1 = 2. * p1 - 2. * p1 * p1
        alpha2 = 2. * p2 - 2. * p2 * p2    
        num = ( p1 -p2)**2 - (n1+n2)*(n1*alpha1 + n2 * alpha2)/(4.*n1*n2*(n1+n2-.1))
        denom = (p1 -p2)**2 + (4.0*n1*n2-n1-n2)*(n1*alpha1 + n2 * alpha2)/(4.*n1*n2*(n1+n2-1.))
        if denom == 0:
            return 0.0
        return num/denom
        
    def getFSTForSets(self):
        fst=[]
        for sfs in self.setMDSFS:
            fst_i=0
            for i in range(self.n1+1):
                for j in range(self.n2+1):
                    fst_i+=sfs[i,j]*self.getFST1Locus(float(i)/self.n1,float(j)/self.n2,self.n1/2.,self.n2/2.)
            fst.append(fst_i/sum(sum(sfs)))
        return fst
    
    def getFSTSimple(self):
        pi1     =   self.getPi(individuals=self.pops[0]) 
        pi2     =   self.getPi(individuals=self.pops[1])
        piw     =   (pi1+pi2)/2

        pib     =   self.getPi()
        
        return 1-piw/pib

    def getFSTSimpleForSets(self):
        self.getPiFromSets()
        self.fst=[]

        for k,pi in enumerate(self.pi):
            piw=(pi[0]+pi[1])/2
            pib=pi[2]
            self.fst.append(1-piw/pib)
        return self.fst

    def runMsSims(self,command,pops=None):
        command +=">data.txt"
        os.system(command)
        self.readMsFile("data.txt")
        self.setPops(pops)
    def runMsSimsC(self,command,pops=None):
        command +=">data.txt"
        os.system(command)
        self.readMsFileCombine("data.txt")
        self.setPops(pops)
        print self.data.shape

    def writeFreqFile(self,file="freq.txt"):
        if not hasattr(self,"mdfreq"):
            self.createMultiDimSFS()
        f=open(file,"w")
        for i in np.transpose(self.mdfreq):
            f.write("%d %d\n"%(i[0],i[1]))
        f.close()

    def writePopSizeFile(self,file="ss.txt"):
        f=open(file,"w")
        for i,p in enumerate(np.bincount(self.popVector)):
            if i>0: f.write("pop%d %d\n"%(i,p))
        f.close()

#---------------------------------------------------------------------------------
    def simOneGen(self,newNHap):
        """old data is safed in self.data. The following approach
        is straightforward, but also inefficient:
            1. simply chose a gamete from prev"""
        if self.nPops   ==  1:
            v           =   rng.randint(0,self.nHap,newNHap)
            newData     =   np.array([ self.data[i] for i in v])
            self.nHap   =   newNHap
            self.data   =   newData
            return

    def simOneGenRec(self,newNHap,recRate,verbose=False):
        """for now, just assume uniform recombination"""
        nRecEvents=rng.poisson(recRate,newNHap)
        if verbose: print nRecEvents
        newData=np.empty((newNHap,self.nSegsites))

        for i,nre in enumerate(nRecEvents):
            if nre==0:
                newData[i]      =   self.data[rng.randint(0,self.nHap)]
                if verbose: print i,"zero"
            else:
                v               =   rng.randint(0,self.nHap,nre+1)
                recPos          =   np.empty(nre+2);recPos[0]=0;recPos[nre+1]=1
                recPos[1:nre+1] =   np.sort(rng.random(nre))

                nChr            =   np.empty(self.nSegsites)
                for j in range(nre+1):
                    pos         =   np.logical_and(self.segSites>recPos[j],
                                    self.segSites<recPos[j+1])
                    nChr[pos]   =   self.data[v[j]][pos]
                newData[i]      =   nChr
                if verbose: print i,v,recPos

        self.nHap   =   newNHap
        self.data   =   newData

    def simOneGenInfRec(self,newNHap):
        newData=np.empty((newNHap,self.nSegsites),dtype="int")
        for i in range(self.nSegsites):
            oldFreq             =   float(sum(self.data[:,i]))/self.nSegsites
            newData[:,i]   =   rng.random(newNHap)<oldFreq 

        self.nHap   =   newNHap
        self.data   =   newData

    def setupSim(theta=0,rho=0):
        self.theta=theta
        self.rho=rho
    def simSurfGen(self,subPopSize=10):
        if not hasattr(self,"theta"):
            self.setupSim()

        newnHap=nHap+subPopSize
        newData=np.empty((nHap+subPopSize,nSegsites))


    def createMDSFSForSets(self,folded=False):
        """folded nyi"""
        self.setMDSFS=[]
        for s in self.sets:
            condFreq=self.mdfreq[:,s]
            sfs=np.zeros((len(self.pops[0])+1,len(self.pops[1])+1),dtype="int")
            for i in np.arange(condFreq.shape[1]):
                sfs[tuple(condFreq[:,i])]+=1
            self.setMDSFS.append(sfs)
        if self.setsIncludeGlobal:
            self.setSFS.insert(0,self.mdsfs)
        self.setMDSFS=np.array(self.setMDSFS)

    def createSFSForSets(self):
        self.setSFS=[]
        for s in self.sets:
            condFreq=self.mdfreq[:,s]
            sfs=[np.zeros(len(self.pops[0])+1),
                 np.zeros(len(self.pops[1])+1),
                 np.zeros(self.nHap)
                ]
            for i in np.arange(condFreq.shape[1]):
                sfs[0][condFreq[0,i]]+=1
                sfs[1][condFreq[1,i]]+=1
                sfs[2][sum(condFreq[:,i])]+=1
            self.setSFS.append(sfs)
        if self.setsIncludeGlobal:
            condFreq=self.mdfreq
            sfs=[np.zeros(len(self.pops[0])+1),
                 np.zeros(len(self.pops[1])+1),
                 np.zeros(nHap)
                ]
            for i in np.arange(condFreq.shape[1]):
                sfs[0][condFreq[0,i]]+=1
                sfs[1][condFreq[1,i]]+=1
                sfs[2][sum(condFreq[:,i])]+=1
            self.setSFS.insert(0,sfs)

    def getPiFromSets(self):
        self.pi=[]

        n1=len(self.pops[0])
        n2=len(self.pops[1])
        ntot=self.nHap

        den1 = n1*(n1-1.0)
        den2 = n2*(n2-1.0)
        dent = ntot*(ntot-1.0)
        
        for k,sfs in enumerate(self.setSFS):
            thetaPi=np.array([0.0,0.0,0.0])
            for i, freq in enumerate(sfs[0]):
                if i==0 or i==len(self.pops[0]):
                    continue
                thetaPi[0]+=2.0*i*(n1-i)*freq

            for i, freq in enumerate(sfs[1]):
                if i==0 or i==len(self.pops[1]):
                    continue
                thetaPi[1]+=2.0*i*(n2-i)*freq

            for i, freq in enumerate(sfs[2]):
                if i==0 or i==self.nHap:
                    continue
                thetaPi[2]+=2.0*i*(ntot-i)*freq
            
            self.pi.append(thetaPi/np.array([den1,den2,dent]))
        return np.array(self.pi)
    
    def makeBlockJKSets(self,windowSize=0.05,offset=0.05,
                        startEndPos=None,folded=False):
        self.setSlidingWindowSets(windowSize,offset,startEndPos,removeEmpty=False)
        self.createMDSFSForSets(folded=folded)
        self.createMultiDimSFS(folded=folded)
        for i in range(len(self.setMDSFS)):
            self.setMDSFS[i]=self.mdsfs-self.setMDSFS[i]
            
        for i in range(len(self.sets)):
            self.sets[i]=np.setdiff1d(np.arange(self.nSegsites),self.sets[i])
    
    def getBlockJKStandardError(self,setf):
        n = float(len(self.sets))
        theta_i=np.array(setf())
        theta_dot=sum(theta_i/n)
        se_jack=np.sqrt((n-1)/n*sum((theta_i-theta_dot)**2))
        return se_jack

    def getBlockJKStandardError2(self,setf,f):
        n=self.nSegsites
        g=float(len(self.sets))
        m=n-np.array([float(len(self.sets[i])) for i in xrange(len(self.sets))])
        h=n/m
        theta_hat_n=f()
        theta_hat_jStar=np.array(setf())

        intSum=sum((1-(m/n))*theta_hat_jStar)
        extSum=sum((1/(h-1))*(h*theta_hat_n-(h-1)*theta_hat_jStar-g*theta_hat_n+intSum)**2)
        return 1/g * extSum


    def doDadi2(self,sfs,freqPath="/data/surfing/dadi_freqs/",folded=False):
        import dadi
        #print sfs
        sfs=dadi.Spectrum(sfs)
        sfs=sfs.fold()
        #print sfs
        def seq(min,max,step=50):
            return (np.arange(0,1+1./step,1./step,dtype="f")*(max-min))+min

        p1=np.exp(seq(np.log(0.1),np.log(50),50))
        p2=np.exp(seq(np.log(0.01),np.log(100),100))
        p3=list(p2)
        if freqPath=="/data/surfing/dadi_freqs/":
            p3.extend([1.1,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09])
        else:
            p3.extend([0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.1])
        ll=np.zeros((51,120))
        for i in range(51):
            for j in range(120):
                model=np.loadtxt("%s/freq_%d_%d.txt"%(freqPath,i,j))
                model=dadi.Spectrum(model)
                if folded:
#                    print "FOLDED",i,j,
                    model.folded=True
                    model.mask=model==np.exp(1)


                ll[i,j]=dadi.Inference.ll_multinom(model,sfs)
#        model=np.array([[dadi.Inference.ll_multinom(dadi.Spectrum(np.loadtxt("/data/surfing/dadi_freqs/freq_%d_%d.txt"%(i,j))),sfs)
#                      for j in range(120)] for i in range(51)])
        maxPos=np.where(ll==np.max(ll))
        print maxPos
        return p3[maxPos[1]]
