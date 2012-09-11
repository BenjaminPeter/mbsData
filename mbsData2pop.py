#!/usr/bin/python

from mbsData import mbsData
import numpy as np
import os
import numpy.random as rng

class mbsData2P(mbsData):
    def __init__(self):
        self.nPops=1
        self.pops       =   None
        self.popVector  =   None

    def setPops(self, pops):
        """pops should be a list of lists of the IDs for each individual in the population"""
        nPops=len(pops)
        self.pops=[np.array(p) for p in pops]
        self.popVector=np.zeros(self.nHap,dtype="i4")
        for i,p in enumerate(pops):
            for pp in p:
                self.popVector[pp]=i+1
        self.n1=len(self.pops[0])
        self.n2=len(self.pops[1])
        #print self.popVector
    
    def createMultiDimSFS(self):
        mdfreq=np.empty((nPops,self.nSegsites),dtype="I")
        for i in range(self.nSegsites):
            snp=self.getDataFromId(id=i)
            for j in range(nPops):
                #print i,j,snp
                mdfreq[j,i]=sum(snp[self.pops[j]])

        dimensions=[len(self.pops[0])+1,len(self.pops[1])+1]
        #print dimensions
        mdsfs=np.zeros(dimensions,dtype="I")

        shared=0
        private1=0
        private2=0
        for i in range(self.nSegsites):
            mdsfs[tuple(mdfreq[:,i])]+=1
            if len(mdsfs.shape)==2:
                if mdfreq[0,i]==self.n1 and mdfreq[1,i]==self.n2:
                    continue
                if mdfreq[0,i]>0 and mdfreq[1,i]>0:
                    shared+=1
                elif mdfreq[0,i]>0:
                    private1+=1
                elif mdfreq[1,i]>0:
                    private2+=1
                else:
                    print "ERROR: neither private nor shared "

        mdsfs[len(self.pops[0]),len(self.pops[1])]=0
        self.shared=shared
        self.private=(private1,private2)
        self.mdfreq,self.mdsfs=mdfreq,mdsfs
        return mdfreq,mdsfs

	def _getSurfStatColumn(self,sfs,l,u,b,shared,verbose):
		"""gets SFS for column i, summed from l to u with constant b, on "shared" shared SNPs"""
		k=0.
		if shared==0:
			return 0
		sfsSum=np.sum(sfs[l:u[0],l:u[1]])
		sfs/=float(sfsSum)
		a1=np.sum(sfs[l:u[0],l:u[1]],axis=1)*range(l,u[0])*b[0]
		a2=np.sum(sfs[l:u[0],l:u[1]],axis=0)*range(l,u[1])*b[1]
		return(a-b)

	
	def _getSurfStatAscertainmentBias(self,sfs,l=1,u=0,b=(1.,1.),verbose=True):
		pass
			
    def _getSurfStatData(self,sfs,l,u,b,shared,verbose,excludeDiag=False):
        k=0.
        for i in range(l,u[0]):
            for j in range(l,u[1]):
                if i==j and excludeDiag:
                    shared -= sfs[i,j]
                else:
                    k+=sfs[i,j]*(b[0]*i-b[1]*j)
                if verbose: print k,i,j,sfs[i,j]
        if shared==0:
            return 0
        return float(k)/shared

    def _getSurfStat(self,l=1,u=0,b=(1.,1.),verbose=True,excludeDiag=False):
        if not hasattr(self,"mdsfs"):
            self.createMultiDimSFS()
        if self.shared==0:
            return 0

        return self._getSurfStatData(self.mdsfs,l,(self.n1+u,self.n2+u),b,
                                     self.shared, verbose,excludeDiag)

    def getPhi1(self,verbose=False,excludeDiag=False):
        return self._getSurfStat(1,1,(1.,1.),verbose,excludeDiag)

    def getPhi2(self,verbose=False):
        n1=self.n1; n2=self.n2
        return self._getSurfStat(1,1,((1.+n2)/(n1*n2),(n1+1.)/(n1*n2)),verbose)

    def getPhi3(self,verbose=False):
        n1=self.n1; n2=self.n2
        return self._getSurfStat(0,0,(1./n1,1./n2),verbose)

    def getPhi4(self,verbose=False):
        n1=self.n1; n2=self.n2
        return self._getSurfStat(1,0,((1.+n2)/(n1*n2),(n1+1.)/(n1*n2)),verbose)
    
    def getPhi1Ex(self,verbose=False):
        return self.getPhi1(verbose,True)

    def _getSurfStatForSets(self,l=1,u=0,b=(1.,1.),
                            verbose=True,excludeDiag=False):
        phi=[]
        for sfs in self.setMDSFS:
            shared=sum(sum(sfs[1:,1:]))
            if shared==0:
                phi.append(0)
            else:
                phi.append(self._getSurfStatData(sfs,l,(self.n1+u,self.n2+u),b,shared,
                                                 verbose,excludeDiag))
        return phi

    def getPhi1ForSets(self,verbose=False,excludeDiag=False):
        return self._getSurfStatForSets(l=1,u=1,b=(1.,1.),verbose=verbose,
                                       excludeDiag=excludeDiag)

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
        num=0.0;denom=0.0
        fst=0.0
        for i in range(self.n1+1):
            for j in range(self.n2+1):
                x=self.getFST1Locus(float(i)/self.n1,float(j)/self.n2,self.n1,self.n2)
                num_c,denom_c=x
                num+=num_c*self.mdsfs[i,j]
                denom+=denom_c*self.mdsfs[i,j]
                #fst += self.mdsfs[i,j]*self.getFST1Locus(float(i)/self.n1,float(j)/self.n2,self.n1/2.,self.n2/2.)
        #return fst/self.nSegsites
        return num/denom

    def getFST1Locus(self,p1,p2,n1,n2):
        alpha1 = 2. * p1 - 2. * p1 * p1
        alpha2 = 2. * p2 - 2. * p2 * p2    
        num = ( p1 -p2)**2 - (n1+n2)*(n1*alpha1 + n2 * alpha2)/(4.*n1*n2*(n1+n2-.1))
        denom = (p1 -p2)**2 + (4.0*n1*n2-n1-n2)*(n1*alpha1 + n2 * alpha2)/(4.*n1*n2*(n1+n2-1.))
#        if denom == 0:
#            return 0.0
        return num,denom


        
    def getFSTForSets(self):
        fst=[]
        for sfs in self.setMDSFS:
            fst_i=0
            num=0.0;denom=0.0
            for i in range(self.n1+1):
                for j in range(self.n2+1):
                    num_c,denom_c=self.getFST1Locus(float(i)/self.n1,float(j)/self.n2,self.n1/2.,self.n2/2.)
                    num=num_c*sfs[i,j]
                    denom=denom_c*sfs[i,j]
            fst_i=num/denom
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


    def createMDSFSForSets(self):
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
    
    def makeBlockJKSets(self,windowSize=0.05,offset=0.05, startEndPos=None):
        self.setSlidingWindowSets(windowSize,offset,startEndPos,removeEmpty=False)
        self.createMDSFSForSets()
        self.createMultiDimSFS()
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

    def readHGDPBlocks(self):
        sets=dict()
        f=open("snp2.txt","r")
        i=0
        while True:
            line=f.readline()
            if line=="":
                break
            chr,pos,id,s=line.split()
            if not int(s) in sets:
                sets[int(s)]=[]
            sets[int(s)].append(i)
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



    
    

    
    


