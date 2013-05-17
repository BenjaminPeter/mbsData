#!/usr/bin/python
import numpy as np
import fileinput
import os

class mbsData:
    def __init__(self,file=None,verbose=False):
        """if file is given, it loads the file assuming mbs data"""
        """if verbose, some stuff is put out"""
        if file != None:
            self.readMbsFile(file)

        self.verbose = verbose

    def _getDefault(self,pos=None,id=None,individuals=None
                       ,singleSite=False,derivedOnly=False):
        """for streamlining all functions. pos, id give position in bp and id of
        SNP, individuals give a subset of individuals to consider"""
        if not hasattr(self,"length"):
            if hasattr(self,"length"):
                self.length=max(self.segSites)+1
            else:
                self.length=1.0
        if pos==None and id == None:
            if hasattr(self,"selPos") and singleSite:
                pos =   self.selPos
                id  =   self.selId
            elif singleSite==False:
                pos =   (0,self.length)
                id  =   (0,self.nSegsites-1)
            else:
                raise NameError("no selected Site defined")
        #only id is set
        elif id==None:
            if isinstance(pos,tuple):
                id   =   (self.getClosestSNPfromPos(pos[0],1),
                          self.getClosestSNPfromPos(pos[1],-1))
            elif type(pos)==np.ndarray or type(pos)==list:
                id  =   np.array([self.getClosestSNPfromPos(i) for i in pos])
                if type(pos)==list:
                    pos =   np.array(pos)
            else:
                id  =   self.getClosestSNPfromPos(pos)
        elif pos==None:
            if isinstance(id,tuple):
                pos  =   (self.getPosFromId(id[0]),
                          self.getPosFromId(id[1]))
            elif type(id)==np.ndarray or type(id)==list:
                pos  =   np.array([self.segSites[i] for i in id])
                if type(id)==list:
                    id  =   np.array(id)
            else:
                pos =   self.getPosFromId(id)


        if individuals==None:
            if derivedOnly == True:
                individuals=self.getIndividualsWithDerivedAllel(id=id)
            else:
                individuals=np.arange(self.nHap)
        individuals=np.array(individuals)
            

        return (pos,id,individuals)
#-----------------------------set predefined intervals----------------------------
    def setSets(self,sets,includeGlobal=True):
        self.sets=sets
        self.setsIncludeGlobal=includeGlobal

    def setSequencedRegionTrack(self,seqFile,format="anna"):
        """this function defines for which region information is available, to
        remove unsequenced region
        file has the following format:
            1. row: header
            2-Nth row: data
            1. col: id (not read)
            2. col type (not read)
            3. col: chromosome (read if chromosome info is available)
            4. col: start of a fragment
            5. col: end of a fragment
        output is saved in self.seqRegion as a tuples of start/endposlist
        """
        """ history: 30.8.2011: created"""
        """ history: 10.10.2011: fixed a bug with readChr"""
        
        #first check if we have chr information
	readChr=False
        if hasattr(self,"chromosome"):
            readChr=True
        
        #read File
        fp = open(seqFile,"r")
        fp.readline() #discard header
        minPos=min(self.segSites)
        maxPos=max(self.segSites)
        self.seqRegion=([],[])
        for l in fp:
            line = l.split()
            if readChr:
                if self.chromosome != line[3]: continue
            il4 = int(line[4])
            il5 = int(line[5])
            if (il4 > minPos and il4 < maxPos) or \
               (il5 > minPos and il5 < maxPos):
                self.seqRegion[0].append(int(line[4]))
                self.seqRegion[1].append(int(line[5]))

    def removeMonomorphicSites(self):
        self.createSFS()
        toKeep = np.where(self.freqs!=0)[0]
        self.segSites=np.array(self.segSites)[toKeep]
	if hasattr(self,"polarizable"): self.polarizable=self.polarizable[toKeep]
        self.data=self.data[:,toKeep]
        self.nSegsites=len(self.segSites)
        self.setSelectedSite(pos=self.selPos)

        del self.SNPdict
        self.createSNPdict()

    def removeSNP(self,id):
        toKeep=np.ones(self.nSegsites,dtype="Bool")
        toKeep[id]=0
        toKeepId = np.where(toKeep)[0]
        self.segSites=np.array(self.segSites)[toKeepId]
	if hasattr(self,"polarizable"): self.polarizable=self.polarizable[toKeepId]
        self.data=self.data[:,toKeepId]
        self.nSegsites=len(self.segSites)
        self.setSelectedSite(pos=self.selPos)


        del self.SNPdict
        self.createSNPdict()
        
    def removeSNPNotSequenced(self):
        """removes the SNP outside of sequenced regions, e.g. because they are
        from a different source. 
        """

        if not hasattr(self,"seqRegion"):
            return

        #get sites that should remain in data
        toKeep=np.zeros(self.nSegsites,dtype="Bool")
	self.segSites=np.array(self.segSites)
        for i,sp in enumerate(self.seqRegion[0]):
            ep = self.seqRegion[1][i]
            kp = np.logical_and(self.segSites>sp ,self.segSites<ep)
            toKeep = np.logical_or(kp,toKeep)

        toKeepId = np.where(toKeep)[0]
        self.segSites=np.array(self.segSites)[toKeepId]
	if hasattr(self,"polarizable"): self.polarizable=self.polarizable[toKeepId]
        self.data=self.data[:,toKeepId]
        self.nSegsites=len(self.segSites)
        self.setSelectedSite(pos=self.selPos)


        del self.SNPdict
        self.createSNPdict()
        
        return toKeep

    def setExpandingWindowSets(self,step=1000,startEndPos=None):
        self.windowSize=step
        if startEndPos == None and hasattr(self,"seqRegion"):
            startEndPos=(self.seqRegion[0][0],self.seqRegion[1][-1])
        elif startEndPos == None:
            startEndPos=(min(self.segSites),max(self.segSites))
        startPos=startEndPos[0]
        endPos=startEndPos[1]
        
        startPoints=np.array([startPos for i in range(startPos,endPos,step)])
        endPoints=np.array([i for i in range(startPos,endPos,step)])
        midPoints=(startPoints+endPoints)/2


        self.sets=[]
        #now get SNP in each window
        for i,sp in enumerate(startPoints):
            ep = endPoints[i]
            self.sets.append(np.where(np.logical_and(self.segSites>sp,self.segSites<ep))[0])

        self.windows=(np.array(startPoints),np.array(endPoints),np.array(midPoints))


        self.setsIncludeGlobal=False
        self.createSFS()
        self.createSFSForSets()

    def setSlidingWindowSets(self,windowSize=10000,offset=5000,startEndPos=None,removeEmpty=False):
        """this function calculates sliding windows and saves them in sets.
        if seqRegion is set, it also calculates for each window how much of it is
        covered by each window
        startEndPos is a tuple(startPos,endPos) which gives the start and end of
        the sliding window. if it set to None, it is inferred
        """
        """history: 30.8.2011: created"""

        self.windowSize=windowSize
        if startEndPos == None and hasattr(self,"seqRegion"):
            startEndPos=(self.seqRegion[0][0],self.seqRegion[1][-1])
        elif startEndPos == None:
            startEndPos=(min(self.segSites),max(self.segSites))
        startPos=startEndPos[0]
        endPos=startEndPos[1]
        
        midPoints=np.arange(startPos+windowSize/2.,endPos-windowSize/2.+0.001,offset)
        startPoints=midPoints-windowSize/2
        endPoints=midPoints+windowSize/2


        self.sets=[]
        #now get SNP in each window
        for i,sp in enumerate(startPoints):
            ep = endPoints[i]
            self.sets.append(np.where(np.logical_and(self.segSites>sp,self.segSites<ep))[0])

        #for each window get the coverage of sequence
        coverage=[]
        if hasattr(self,"seqRegion"):
            nSeqRegion=len(self.seqRegion[0])
            c=0 #current seqRegion
            for i,startWin in enumerate(startPoints):
                endWin = endPoints[i]
                cv=0
                for j,startSeq in enumerate(self.seqRegion[0]):
                    endSeq=self.seqRegion[1][j]


                    #sequencing ends  before window or starts after window
                    if endSeq < startWin or startSeq > endWin:
                        #print "result: no overlap",cv
                        continue

                    #if sequence is fully contained in window, add full seq size
                    elif startWin <= startSeq and endSeq <= endWin:
                        cv+=endSeq-startSeq
			print "[%i - %i] vs [%i -%i]" %\
                        (startWin,endWin,startSeq,endSeq),
                        print "result: fully contained",cv
                    elif startSeq <= startWin and endWin <= endSeq:
                        cv+=endWin -startWin
                    #else, right overlap
                    elif startSeq <= startWin and endSeq <= endWin:
                        cv+=endSeq-startWin
			print "[%i - %i] vs [%i -%i]" %\
                        (startWin,endWin,startSeq,endSeq),
                        print "result: right overlap",cv
                    #else, left overlap
                    elif startWin < startSeq and endWin < endSeq:
                        cv+=endWin-startSeq
                        print "[%i - %i] vs [%i -%i]" %\
                        (startWin,endWin,startSeq,endSeq),
                        print "result: left overlap",cv
                    else:
                        pass#print "WAAAAAAAAAAAAAAAAAH", cv
                coverage.append(cv)
        else: coverage = np.array([windowSize for i in midPoints])

        #finally,remove windows with 0 coverage or 0 snp
        toRemove=[]
        for i,s in enumerate(self.sets):
            cv = coverage[i]
            if cv == 0 or len(s) ==0:
                toRemove.append(i)

        midPoints = midPoints.tolist()
        startPoints = startPoints.tolist()
        endPoints = endPoints.tolist()
        if removeEmpty:
            for i in toRemove[::-1]:
                self.sets.pop(i)
                coverage.pop(i)
                midPoints.pop(i)
                startPoints.pop(i)
                endPoints.pop(i)

        self.coverage=np.array(coverage)
        self.windows=(np.array(startPoints),np.array(endPoints),np.array(midPoints))


        self.setsIncludeGlobal=False
        self.createSFS()
        self.createSFSForSets()

#-----------------------------make subset----------------------------
    def getSubset(self,pos=None,id=None,individuals=None,setNewSelSite=False):
        """gets a subset of the data as a new mbsData object. if setNewSelSite
        is true, if the selected Site is not in the subset, the next best site
        is put as selected site"""
        pos,id,individuals=self._getDefault(pos,id,individuals)

        if isinstance(id,tuple):
            id=np.arange(id[0],id[1]+1)

        subset = mbsData()
        #stuff that all data sets should have
        subset.data=self.data[individuals,:][:,id]
        subset.individuals = np.arange(len(individuals))
        subset.length=max(self.segSites)
        subset.currDataSet=1
        subset.nDataSets=1
        subset.nHap,subset.nSegsites=subset.data.shape
        subset.segSites=np.array(self.segSites)[id]
        subset.type=self.type
        subset.verbose=self.verbose

        #selected Site
        if sum(self.selId==id)==1:
            subset.selId=np.where(self.selId==id)[0]
            subset.selPos=self.selPos
            subset.selSiteData=np.array(self.selSiteData)[individuals]
        elif setNewSelSite:
            subset.selId=subset.getClosestSNPfromPos(self.selPos)
            subset.selPos=subset.getPosFromId(subset.selId)
            subset.selSiteData=subset.getDataFromId(subset.selId)

            
        #stuff specific to beagle file
        subset.major=self.major[id]
        subset.minor=self.minor[id]
        subset.individualIDs=np.array(self.individualIDs)[individuals]
                               
        subset.han=np.where([s[:2]=="NA" for s in subset.individualIDs])[0]
        subset.tib = tib2=np.where([s[:2]=="DR" or s[:2]=="NQ" for s in
                                                            subset.individualIDs])[0]

        return(subset)

    def getTib(self):
        return self.getSubset(individuals=self.tib,setNewSelSite=True)

    def getHan(self):
        return self.getSubset(individuals=self.han,setNewSelSite=True)
#-----------------------------set selSite----------------------------
    def setSelectedSite(self,pos=None,id=None):
        pos,id,_    = self._getDefault(pos,id,0,singleSite=True)
        self.selId=id
        self.selPos=pos
        self.selSiteData=self.getDataFromId(id)
#-----------------------------getters----------------------------
    def getData(self, pos=None, id=None):
        """get the SNP data from either a position(default) or ID"""
        pos,id,_    = self._getDefault(pos,id,0,singleSite=True)
        if pos!= None:
            return self.getDataFromPos(pos)
        else:
            return self.getDataFromId(id)

    def getDataFromPos(self,pos):
        """gets the SNP from a position. Throws an error if the position is not segregating"""
        if not hasattr(self,'SNPdict'):
            self.createSNPdict()
        try:
            id = self.SNPdict[pos]
        except:
            raise ValueError('not a segregating site')
        return self.data[:,id]

    def getDataFromId(self,id,individuals=None):
        """gets the SNP with a given id"""
        if individuals==None:
            individuals=self.getIndividuals()
        return self.data[individuals,id]

    def createSNPdict(self):
        """creates a dictionary that maps positions to the corresponding ids"""
        self.SNPdict=dict()
        for i in range(len(self.segSites)):
            self.SNPdict[self.segSites[i]] = i

    def getIDfromPos(self,pos):
        if not hasattr(self,"SNPdict"):
            print "made dict"
            self.createSNPdict()
        return self.SNPdict[pos]

    def getPosFromId(self,id):
        return self.segSites[id]

    def getClosestSNPfromPos(self,pos,direction=0):
        """gives a SNP close to the position. If direction==None or direction
        ==0, the closest SNP
        is returned. is direction==1, the next larger SNP is returned, if
        direction== =1, the next SNP with lower position is returned """

        #lookup SNP, and get from dict if it exists: O(1)?
        if not hasattr(self,'SNPdict'):
            self.createSNPdict()
        try:
            s = self.SNPdict[pos]
            return s
        except:
            pass

        #do binary search: O(lg n)
        segs= self.segSites
        l = len(segs)
        while l>1:
            pivot = l/2
            if pos> segs[pivot]: 
                segs = segs[pivot:]
            else:
                segs = segs[:pivot]
            l= len(segs)
        lower = segs[0]
        lPos = self.SNPdict[lower]
        if lPos ==0:
            if direction==-1:
                raise ValueError('no smaller SNP than this')
            elif pos < self.segSites[0]:
                return 0 
            else:
                return 1


        uPos = lPos+1
        if uPos ==len(self.segSites):
            if direction==1:
                raise ValueError('no SNP genotyped at larger pos')
            else:
                return lPos
        upper = self.segSites[uPos]

        #print lPos, uPos, "diff: ", pos-lower, upper-pos
        if abs(pos - lower) > abs(upper-pos) and direction != -1:
            return uPos
        elif direction != 1:
            return lPos
        return uPos

    def getHaplotype(self,pos=None,id=None,individuals=None):
        """gets a given haplotype, with optional starting and ending position, for the individuals given in individual"""
        pos,id,individuals=self._getDefault(pos,id,individuals,singleSite=False)

        if type(id)==np.ndarray:
            return np.array(self.segSites)[id],self.data[individuals,:][:,id]


        if id[0]    ==  id[1]: #only one SNP
            return self.segSites[id[0]],self.data[individuals,id[0]]

        return self.segSites[id[0]:id[1]+1],self.data[individuals,id[0]:id[1]+1]
    
    def getIndividualsWithDerivedAllel(self,pos=None,id=None,individuals=None):
        """gets all the individuals with the derived allel at the given snp. If id==None, the selected site is assumed"""
        pos,id,individuals=self._getDefault(pos,id,individuals,singleSite=True)
        return np.where(self.getDataFromId(id)==1)[0]

    def getIndividualsWithAncestralAllel(self,pos=None,id=None,individuals=None):
        """gets all the individuals with the ancestral allel at the given snp. If id==None, the selected site is assumed"""
        pos,id,individuals=self._getDefault(pos,id,individuals,singleSite=True)
        return np.where(self.getDataFromId(id)==0)[0]

    def getIndividuals(self):
        """gets all the individuals in the sample"""
        return np.arange(self.nHap)

#-----------------------------mason----------------------------
    def getNSL(self):
        sl_max=200
        def pwcomparisons(mutants,sn):
            monster=np.ones((self.nSegsites,self.nHap,self.nHap)) #a zeros or ones? a key question.
            for s in range(self.nSegsites):
                for m in mutants[s]:
                    monster[s,m,:]=1-monster[s,m,:]
                    monster[s,:,m]=1-monster[s,:,m]
            return monster
        sn=self.nSegsites
        n=self.nHap

        d2=zip(*self.data)
        m2=map(lambda x:[y for y in range(n) if x[y]],d2)
        monster=pwcomparisons(m2,self.nSegsites)
        forward=np.zeros((sn,n,n))
        backward=np.zeros((sn,n,n))
        forward[0]*=monster[0]
        backward[-1]*=monster[-1]
        for s in range(1,sn):
            forward[s]=forward[s-1]+1
            forward[s]*=monster[s]
            backward[-1-s]=backward[-s]+1
            backward[-1-s]*=monster[-1-s]

        forward=np.roll(forward,1,axis=0)  
        sl=forward+backward+1
        sl=np.minimum(sl,sl_max)
        sl*=np.tile([np.tri(n,k=-1)],(sn,1,1))

        anc=np.ones((sn,n,n))
        der=np.zeros((sn,n,n))
        for s in range(sn):
            for m in m2[s]:
                anc[s,m,:]=0
                anc[s,:,m]=0
                der[s,m,:]+=1
                der[s,:,m]+=1
        der=np.maximum(der-1,0)
        anc*=sl
        der*=sl

        daf=np.array(map(len,m2))
        daf[daf==1]=-1
        daf[daf==n-1]=-1 
        anc=anc.sum(axis=1).sum(axis=1)/((n-daf)*(n-daf-1)/2)
        der=der.sum(axis=1).sum(axis=1)/(daf*(daf-1)/2)
        self.nsl=daf,anc,der
        return daf[self.selId],anc[self.selId],der[self.selId]
        #return daf, anc, der
    def getNSLExternal(self):
        self.writeToMsFile("temp.ms")
        s="./nsl -msfile temp.ms >temp.txt 2> /dev/null"
        #print s
        c=0
        while c<100:
            c +=1
            os.system(s)
            fp = open("temp.txt","r")
            fp.readline()

            while 1:
                line = fp.readline()
                line = line.split()
                if len(line)<3:
                    break
                if float(line[2])==0.0:
                    #print "Null"
                    break

                if int(line[1]) == self.selPos:
                    return float(line[2])
                if line=="\n":
                    break
        return -100000
        

    def getNSLforSetsAvg(self,requirePolarizedData=False):
        if not hasattr(self,"nsl"):
            self.getNSL()
        nsl_anc=np.ma.array(self.nsl[1])
        nsl_des=np.ma.array(self.nsl[2])
        nsl=nsl_anc/nsl_des
        nsl[np.logical_or(nsl==np.inf,-nsl==np.inf)]=np.ma.masked
        if requirePolarizedData:
            nsl[np.where(np.logical_not(self.polarizable))[0]] = np.ma.masked

        nslSets=np.zeros(len(self.sets))
        for i,set in enumerate(self.sets):
            nslSets[i]=np.mean(nsl[set])
        return nslSets



#-----------------------------IHS/EHH----------------------------
    def getEHHforSetsAvg(self,x):
        """"""
        #first, calculate all EHH values
        ehh=np.empty(self.nSegsites)
        for i,ss in enumerate(self.segSites):
            ehh[i]=self.getEHH(x,id=i,derivedOnly=True)

        #then get avg for sets
        ehhSets=np.empty(len(self.sets))
        for i,set in enumerate(self.sets):
            ehhSets[i]=np.mean(ehh[set])
        self.ehh=ehhSets
        return self.ehh

    def getIHSforSetsAvg(self,requirePolarizedData=False):
        ihs=np.empty(self.nSegsites)
        for i,ss in enumerate(self.segSites):
            ihs[i]=self.getIHS(id=i)[0]
        ihs=np.ma.array(ihs)
        ihs[np.logical_or(ihs==np.inf,-ihs==np.inf)]=np.ma.masked
        if requirePolarizedData:
            ihs[np.logical_not(self.polarizable)] = np.ma.masked


        ihsSets=np.empty(len(self.sets))
        for i,set in enumerate(self.sets):
            
            ihsSets[i]=np.mean(ihs[set])
        self.ihs=ihsSets
        return self.ihs

            
    def _getUniqueHaplotypes(self,data=None,pos=None,id=None,individuals=None):
        """for the data given in data, computes the unique haplotype, returns a dict"""
        if data==None:
            pos,id,individuals=self._getDefault(pos,id,individuals)
            if isinstance(id,tuple):
                id=np.arange(id[0],id[1]+1)

            data=self.data[individuals,:][:,id]
        nHap=data.shape[0]
        unique=dict()
        #if we only have 1 segregating site, return allele frequencies
        if len(data.shape)==1:
            unique[0]=data.tolist().count(0)
            unique[1]=data.tolist().count(1)
            return unique

        #else, change to string... 
        for i in range(nHap):
            s = tuple(j for j in data[i])
            if unique.has_key(s):
                unique[s]+=1
            else:
                unique[s] = 1

        return unique

    def _getEHHsingle(self,x,id,individuals,onesided=0):
        """gets EHH for a single x. Use getEHH instead."""

        pos = self.segSites[id]
        """call get EHH instead"""
        if individuals.shape == (0,):
            return 0 
        if onesided > 0:
            apos,hapData=self.getHaplotype((pos,pos+x),individuals=individuals)
        elif onesided < 0:
            apos,hapData=self.getHaplotype((pos-x,pos),individuals=individuals)
        else:
            apos,hapData=self.getHaplotype((pos-x,pos+x),individuals=individuals)
        #print pos,hapData
        if hapData == []:
            raise ValueError("empty data; ind="+str(individuals))
#        print apos
        return self._getEHHgeneral(hapData)

    def getEHHfromPos(self,x,pos,individuals=None,onesided=0):
        """gets the EHH for a SNP at a given position"""
        id = self.getClosestSNPfromPos(pos)
        if individuals==None:
            ind=self.getIndividualsWithDerivedAllel(id=id)
        else:
            ind=individuals
        return self.getEHH(x,onesided=onesided,individuals=ind,id=id)

    def getEHH(self,x,pos=None,id=None,individuals=None,onesided=0,verbose=None,derivedOnly=False):
        """calculates EHH in the range SelectedSite +/- x sites
                - onesided: if negative, only upstream allels are considered, if
                positive, only downstream allels are considered, otherwise both
                are considered(default)
                - id: the id of the SNP, if None, the selected SNP is assumed
                - individuals: if only a subset of individuals should be used:
                    default=all individuals
                    - verbose: adds some output"""
        pos,id,individuals=self._getDefault(pos,id,individuals,singleSite=True,derivedOnly=derivedOnly)

        if verbose==None:
            verbose = self.verbose
        ehh=[]
        if isinstance(x,tuple) or isinstance(x,np.ndarray):
            for xx in x:
                ehh.append(self._getEHHsingle(xx,id,individuals,onesided))
            return ehh
        else:
            return self._getEHHsingle(x,individuals=individuals,id=id,onesided=onesided)

    def _getEHHgeneral(self,hapData):
        """don't call this. call getEHH instead
            calculates EHH on the haplotypes given in hapData:
            hapData is expected to be an np.array of ints/bools where each row is a haplotype and each column is an individual
        """
        uniqueHaplotypes=self._getUniqueHaplotypes(hapData)
        #print uniqueHaplotypes.values()
        cnt = uniqueHaplotypes.values()
        
        #get Denominator 
        denom = sum(cnt) * (sum(cnt)-1) 
        if denom == 0:
            return 0
            
        #print "denom: ",denom
        try:
            cnt = [cnt.count(i+1) for i in range(max(cnt))]
        except:
            raise ValueError(str(uniqueHaplotypes))


        nominator=0
        #print "found %i haplotypes at freq %i"% (cnt[0],1)
        for i in range(1,len(cnt)):
            nominator+=cnt[i]*i*(i+1)
            #print "found %i haplotypes at freq %i"% (cnt[i],i)

        return float(nominator)/denom

    def getIHSfromPos(self,pos,threshold=0.05,verbose=False):
        """gets IHS from the closest SNP from a given position"""
        id=self.getClosestSNPfromPos(pos)
        return self.getIHS(id,threshold,verbose)

    def _integrateEHH(self,id,individuals,dir,threshold,maxDist=None,
                         interpolation=False, verbose=False):
        """calculates the integrated EHH
            - id = id of the SNP to start at
            - dir = direction, if negative, upstream, if positive, downstream
            - ind = the individuals to consider
            - threshold = threshold where to stop
            - interpolation = if true, the interpolation by Voight et al. is used, if false, a step function more appropriate for sequence data is used
            -maxDist: the maximal distance to integrate to
        """
        ind=individuals 

        pos=self.segSites[id]
        #check if there are more than 2 individuals, otherwise return 0
        if len(ind)<2:
            return 0,-1#pos
        #check if we have the topmost snp
        if id==self.nSegsites-1 and dir>0: 
            return 0,self.segSites[id]#pos
        if id==0 and dir<0:
            return 0,self.segSites[id]

        iEHH=0
        if dir >0:
            dir =1
        else:
            dir = -1

        #if no maxDist is set, set it to +/- infinity
        if maxDist == None:
            maxDist=np.inf
        maxDist = maxDist * dir
        endId               =   self.getClosestSNPfromPos(pos+maxDist,direction=-dir)
        curId,curPos,curEHH =   id,self.segSites[id],1
        newId,newPos           =   id+dir,self.segSites[id+dir]


        while (newId>0 and dir<0) or \
              (newId <self.nSegsites-1 and dir>0):
            delta = (newPos - curPos)*dir
            newEHH = self.getEHH(abs(pos-newPos),id=id,onesided=dir,individuals=ind)

            if interpolation:
                k = (newEHH+curEHH-.1)*.5*delta
            else:
                k = curEHH*delta #works because curEHH is always >= newEHH

            if verbose:
                print "[pos[%f-%f];id[%i/%i]step: %f/ ehh[i]%f/added: %f/total: %f]"  %(newPos,curPos,newId,curId,delta,newEHH,k,iEHH)
            if newEHH < threshold:
                if interpolation:
                    iEHH        +=  (newPos-curPos) * (1 - (0.05-newEHH)/                       \
                                    (curEHH-newEHH))*.5*(curEHH-0.05)
                    if verbose: print (newPos-curPos)*(1-(0.05-newEHH)/(curEHH-newEHH))*.5*(curEHH-0.05)
                    print iEHH
                    print "end"
                else:
                    iEHH+=curEHH*delta
                break
            iEHH+=k
            #print "[%f][%f-%f][%f/%f]" %(pos,curPos,newPos,curEHH,newEHH)
            curPos,curId,curEHH         =       newPos,newId,newEHH
            newPos,newId                =       self.segSites[newId+dir],newId+dir


            #do some adjustment for distance between last snp and end of ehh
            if curId == endId:
                if curId == 0 or curId == self.nSegsites -1:
                    break
                else:
                    #            xmax    -   
                    dist        =   dir * maxDist + dir*(pos-curPos)
                    iEHH += dist*newEHH
                    newPos=(pos+maxDist)
                    break

 
        return iEHH,newPos
            
    def getIHSnoAncDer(self,pos=None,id=None,individuals=None,
                       threshold=.05,verbose=False,
                       interpolation=False,maxDist=None):
        """gets the integrated IHS for all the individuals, ignoring
        ancestral/derived alleles"""
        pos,id,individuals      =       self._getDefault(pos,id,individuals,singleSite=True)
        ihs,dmin                =       self._integrateEHH(id,individuals,-1,threshold,
                                                           maxDist,interpolation,verbose)
        k,dmax                  =       self._integrateEHH(id,individuals,1,threshold,
                                                           maxDist,interpolation,verbose)
        ihs                     +=      k
        return ihs,(dmin,dmax)

    def getIHS(self,pos=None,
               id=None,individuals=None,threshold=0.05,maxDist=None, verbose=False,interpolation=False,noAncDer=False):
        """calculates unstandardized IHS statistic, see Voight et al 2006
            - id = id of the SNP to start at
            - dir = direction, if negative, upstream, if positive, downstream
            - ind = the individuals to consider
            - threshold = threshold where to stop
            - interpolation = if true, the interpolation by Voight et al. is used, if false, a step function more appropriate for sequence data is used
            - if it should just be run once
        """

        if noAncDer:
            return self.getIHSnoAncDer(pos,id,individuals,threshold,verbose,interpolation)
        pos,id,individuals=self._getDefault(pos,id,individuals,singleSite=True)

        iAnc = self.getIndividualsWithAncestralAllel(id=id)
        iDer = self.getIndividualsWithDerivedAllel(id=id)

        ihs_derived,dmin = self._integrateEHH(id,iDer,-1,threshold,maxDist,interpolation,verbose)
        k,dmax= self._integrateEHH(id,iDer,1,threshold,maxDist,interpolation,verbose)
        ihs_derived+=k

        ihs_ancestral,amin= self._integrateEHH(id,iAnc,-1,threshold,maxDist,interpolation,verbose)
        k,amax= self._integrateEHH(id,iAnc,1,threshold,maxDist,interpolation,verbose)
        ihs_ancestral+=k

        try:
            ihs = np.log(ihs_ancestral/ihs_derived)
        except ZeroDivisionError:
            ihs = np.inf

        return ihs,ihs_ancestral,ihs_derived,(amin,amax),(dmin,dmax)

    def getIHSExternal(self):
        np.savetxt("segSites.txt",self.segSites,fmt="%f")
        np.savetxt("dump.txt",self.data,fmt="%i")

        s="./ihs3 segSites.txt dump.txt %i > temp.txt" %self.selId
        #print s
        os.system(s)
	try:
        k=np.loadtxt("temp.txt")
	except IOError:
        k=(np.nan,np.nan,np.nan)
    if k[0] ==-1e13:
        k=(np.nan,np.nan,np.nan)
    return (k[2],k[1],k[0])

    def getXPEHH(self,x,pos=None,id=None,individuals=None,interpolation=False,verbose=False):
        pos,id,individuals = self._getDefault(pos,id,individuals,singleSite=True)

        iPop0 = np.where(self.pops==0)[0]
        iPop1 = np.where(self.pops==1)[0]

        iEHH0,_       =   self._integrateEHH(id,iPop0,-1,0,x,interpolation,verbose)
        k,_           =   self._integrateEHH(id,iPop0,1,0,x,interpolation,verbose)
        iEHH0+=k


        iEHH1,_       =   self._integrateEHH(id,iPop1,-1,0,x,interpolation,verbose)
        k,_           =   self._integrateEHH(id,iPop1,1,0,x,interpolation,verbose)
        iEHH0+=k

        return(np.log(iEHH0/iEHH1),iEHH0,iEHH1)
        
#-----------------------------Heterozygosity----------------------------
    def getHFromPos(self,interval=None):
        if interval ==None:
            interval=(1,self.length)
        if type(interval) != tuple:
            return self.getHFromSNP(self.getClosestSNPfromPos(interval))
        interval=(self.getClosestSNPfromPos(interval[0],1),
                      self.getClosestSNPfromPos(interval[1],-1))
        return self.getH(interval)

    def getHFromSets(self):
        if not hasattr(self,"pi"):
            self.getPiFromSets()

        if not hasattr(self,"S"):
            self.getSFromSets()

        return self.pi/self.S

        

    def getH(self,pos=None,id=None,individuals=None):
        pos,id,individuals=self._getDefault(pos,id,individuals,singleSite=False)
        interval=id

        if isinstance(interval,tuple):
            if interval[1]-interval[0] == 0:
                return self.getHFromSNP(interval[0])
            interval=np.arange(interval[0],interval[1]+1)

        
        H=0


        for i in np.array(interval):
            H+=self.getHFromSNP(i,individuals)
        H/=len(interval)
        return H
        
    def getHFromSNP(self,id,individuals=None):
        if individuals==None:
            individuals = self.getIndividuals()
        SNP = self.getDataFromId(id,individuals)
        p =  float(sum(SNP))
        N = float(self.nHap)
        H = 2*p*(N-p)/(N*(N-1))
        return(H)
#-----------------------------Pi----------------------------
    def getPiForTwoHaplotypesFromId(self,ind1,ind2,interval=None):
        if interval==None:
            interval=(0,self.length)
        _,h1=self.getHaplotype(id=interval,individuals=ind1)
        _,h2=self.getHaplotype(id=interval,individuals=ind2)
        pwd=0
        for i in range(len(h1)):
            if h1[i]!=h2[i]:
                pwd+=1
        return float(pwd)

    def getPiSlow(self,individuals=None,interval=None):
        """calculates pi by pairwise comparisons"""
        if individuals ==None:
            individuals=self.getIndividuals()
        if interval==None:
            interval=(0,self.length)
        interval=(self.getClosestSNPfromPos(interval[0],1),
                      self.getClosestSNPfromPos(interval[1],-1))
        if self.verbose: print interval
        pi=0.0
        l = len(individuals)
        for i in range(l):
            for j in range(i+1,l):
                pi+=self.getPiForTwoHaplotypesFromId(individuals[i],individuals[j],interval)
                if self.verbose: print i,j,pi
        print pi
        return pi*2/l/(l-1)

    def getPiFromSets(self):
        if self.setsIncludeGlobal:
            self.pi=np.empty(len(self.sets)+1)
        else:
            self.pi=np.empty(len(self.sets))

        n=self.nHap
        denominator = n*(n-1.0)

        if self.setsIncludeGlobal:
            self.pi=np.empty(len(self.sets)+1)
        else:
            self.pi=np.empty(len(self.sets))

        for k,sfs in enumerate(self.setSFS):
            thetaPi=0.0
            for i,freq in enumerate(sfs):
                if i==0 or i==self.nHap:
                    continue
                thetaPi+= 2.0*i*(n-i)*freq
            self.pi[k] = thetaPi/denominator
        return self.pi


    def getPi(self,pos=None,id=None,individuals=None):
        """calculates pi from the SFS"""
        pos,id,individuals=self._getDefault(pos,id,individuals)
        interval=id

        if isinstance(interval,tuple):
            if interval[1]-interval[0] == 0:
                return self.getPi(id=interval[0])
            interval=np.arange(interval[0],interval[1]+1)

        pi=0.0
        l = len(individuals)
        nnm1 = l/(l-1.0)

        if type(interval)==np.ndarray:
            arr=interval
        else:
            arr=np.array((interval,))
        for s in arr:
            p = float(sum(self.getDataFromId(s,individuals)))/l
            pi += 2.0*(p*(1.0-p))*nnm1

        return pi

#----------------------------Singletons/sfs--------------------------
    def createSFSForSets(self):
        self.setSFS=[]
        for s in self.sets:
            condFreq=self.freqs[s]
            sfs=np.zeros(self.nHap+1)
            for i in condFreq:
                sfs[i]+=1
            self.setSFS.append(sfs)
        if self.setsIncludeGlobal:
            self.setSFS.insert(0,self.sfs)
        self.setSFS=np.array(self.setSFS)

    def createSFS(self,pos=None,id=None,individuals=None):
        if individuals==None and id==None and pos ==None:
            glob=True
        else:
            glob=False

        pos,id,individuals=self._getDefault(pos,id,individuals)
        interval=id

        if isinstance(interval,tuple):
            interval=np.arange(interval[0],interval[1]+1)

        nHap=len(individuals)
        sfs=np.zeros(nHap+1)
        freqs=np.zeros(self.nSegsites,dtype=float)
        for i in np.array(interval):
            snp=self.getDataFromId(id=i,individuals=individuals)
            p=sum(snp)
            sfs[p]+=1
            freqs[i]=float(p)#/nHap

        if glob:
            self.sfs=sfs
            self.freqs=freqs
        sfs=sfs[1:nHap]
        
        return sfs,freqs
    
    def getNSingletonsFromSets(self):
        return self.setSFS[:,1]

    def getNSingletons(self,pos=None,id=None,individuals=None,folded=True):
        pos,id,individuals=self._getDefault(pos,id,individuals)
        interval=pos

        sfs,_= self.createSFS(individuals=individuals,pos=interval) 
        if folded:
            return sfs[0]+sfs[self.nHap-2]
        else:
            return sfs[0]

    def removeSingletons(self,folded=True):
        self.removeAllelsWithFrequency(freq=1,folded=folded)
    def removeDoubletons(self,folded=True):
        self.removeAllelsWithFrequency(freq=2,folded=folded)
    def removeAllelsWithFrequency(self,freq=1,folded=True):
        if not hasattr(self,"freqs"):
            self.sfs,self.freqs=self.createSFS()
        toKeep=np.where(np.logical_and(self.freqs!=freq ,
                        np.logical_or(folded, self.nHap-self.freqs!=freq)))[0]
        self.data=self.data[:,toKeep]
        self.segSites=np.array(self.segSites)[toKeep]
        self.nSegsites=np.array(self.data.shape)[1]
        self.selId=sum(toKeep<self.selId)

        if hasattr(self,"SNPdict"):
            del self.SNPdict 
            
        del self.sfs, self.freqs
       
        #check if selected site also has to be removed
#        if sum(self.selSiteData)==freq or sum(self.selSiteData) == self.nHap-freq:
#            del self.selSiteData, self.selId, self.selPos

        
#-----------------------------S/Tajima's D/Fay+Wu H----------------------------
    def getFayWuHFromSets(self,requirePolarizedData=False):
        if requirePolarizedData:
            #copy sets
            tempSets=[s.copy() for s in self.sets]

            #polarizableId=np.where(self.polarizable)[0]
            for i,s in enumerate(self.sets):
                self.sets[i]=np.ma.array(s)
                for j,snp in enumerate(s):
                    if not self.polarizable[snp]:
                        self.sets[i][j]=np.ma.masked
            self.sets=[x.compressed() for x in self.sets]
            if hasattr(self,"pi"):
                del self.pi
            if hasattr(self,"setSFS"):
                del self.setSFS
            if hasattr(self,"freqs"):
                del self.freqs
            if hasattr(self,"sfs"):
                del self.sfs
            self.createSFS()
            self.createSFSForSets()
            result = self.getFayWuHFromSets()

            #no undo all the stuff
            del self.pi
            del self.setSFS
            del self.freqs
            del self.sfs
            self.sets=tempSets
            self.createSFS()
            self.createSFSForSets()
            return result
        else:
            if not hasattr(self,"pi"):
                self.getPiFromSets()
            thetaH=0.0
            n=self.nHap
            denominator = n*(n-1.0)

            if self.setsIncludeGlobal:
                self.fwh=np.empty(len(self.sets)+1)
            else:
                self.fwh=np.empty(len(self.sets))

            for k,sfs in enumerate(self.setSFS):
                thetaH=0.0
                for i,freq in enumerate(sfs):
                    if i==0 or i==self.nHap:
                        continue
                    thetaH+= 2*i*i*freq

                self.fwh[k] = self.pi[k] - thetaH/denominator
            return self.fwh
    
    def getFayWuH(self,pos=None,id=None,individuals=None,requirePolarizedData=False):
        pos,id,individuals=self._getDefault(pos,id,individuals)
        interval=id
    
        if isinstance(interval,tuple):
            if interval[1]-interval[0] == 0:
                return self.getHFromSNP(interval[0])
            interval=np.arange(interval[0],interval[1]+1)
        pi=self.getPi(individuals=individuals,id=interval)
        thetaH=0.0
        l=len(individuals)
        nnm1 = 1/l/(l-1.0)

        if type(interval)==np.ndarray:
            arr=interval
        else:
            arr=np.array((interval,))
        for s in arr:
            #check if we want to use SNP
            if (requirePolarizedData and self.polarizable[s]) or \
               not requirePolarizedData:
                p = float(sum(self.getDataFromId(s,individuals)))/l
                thetaH += 2.0*p*p*nnm1
        return pi-thetaH

    def getSFromSets(self):
        if self.setsIncludeGlobal:
            self.S=np.empty(len(self.sets)+1,dtype="I")
        else:
            self.S=np.empty(len(self.sets),dtype="I")
        
        for k,sfs in enumerate(self.setSFS):
            self.S[k] = sum(self.setSFS[k][1:]) 
        return self.S

    def getS(self,pos=None,id=None,individuals=None):
        pos,id,individuals=self._getDefault(pos,id,individuals)
        interval=id
                    
        #this is the number of segregating sites in the entire sample, now
        #filter all sites that are not segregating in the sumbsample
        if isinstance(interval,tuple):
            if interval[1]-interval[0] == 0:
                return self.getHFromSNP(interval[0])
            interval=np.arange(interval[0],interval[1]+1)

        nSegSites=len(interval)
        #for i in range(intervalId[0],intervalId[1]+1):
        for i in np.array(interval):
            snp=self.getDataFromId(i,individuals)
            if np.all(snp) or not np.any(snp):
               nSegSites -= 1

        return nSegSites

    def getTajimasDfromSets(self):
        """ """
        n=self.nHap

        if not hasattr(self,"pi"):
            self.getPiFromSets()

        if not hasattr(self,"S"):
            self.getSFromSets()

        if self.setsIncludeGlobal:
            self.td=np.empty(len(self.sets)+1)
        else:
            self.td=np.empty(len(self.sets))

        #coefficient a1
        a1=0.0
        for i in range(1,n):
            a1+=1./i
        #coefficient a2
        a2=0.
        for i in range(1,n):
            a2+=1./(i*i)

        #more coeffs
        e1 = 1./a1 * ((n+1.)/(3.*(n-1.))-1/a1)
        e2 = 1./(a1*a1+a2) * (((2.*(n*n+n-3.))/(9.*n*(n-1))) - ((n+2.)/(n*a1))+a2/(a1*a1))
        #done with coeff
        
        for k,sfs in enumerate(self.setSFS):
            tdRaw=self.pi[k] - self.S[k]/a1
            vartD=e1*self.S[k]+e2*self.S[k]*(self.S[k]-1.0)
            self.td[k]=tdRaw/np.sqrt(vartD)

        return self.td

    def getTajimasD(self,pos=None,id=None,individuals=None):
        """calculates Tajimas D after the formula in Wakeley 2008 p114/115
           checked against the stats program from ms
        """
        pos,id,individuals=self._getDefault(pos,id,individuals)
        interval=id

        a1=0.0
        n=len(individuals)
        for i in range(1,n):
            a1+=1./i
        #the raw tajimasD
        pi=self.getPi(individuals=individuals,id=interval)
        S=self.getS(individuals=individuals,id=interval)
        tajimasDraw=pi-S/a1

        #first, calculate coefficients a2,e1 and e2
        a2=0.
        for i in range(1,n):
            a2+=1./(i*i)
        e1 = 1./a1 * ((n+1.)/(3.*(n-1.))-1/a1)
        e2 = 1./(a1*a1+a2) * (((2.*(n*n+n-3.))/(9.*n*(n-1))) - ((n+2.)/(n*a1))+a2/(a1*a1))
        vartD=e1*S+e2*S*(S-1.)
        return tajimasDraw/np.sqrt(vartD),S,pi

    def getTD(self,pos=None,id=None,individuals=None):
        t,_,_ = self.getTajimasD(pos,id,individuals)
        return t
#-----------------------------NJ----------------------------
    def _getDistance(self,i1,i2):
        return(sum(self.data[i1]!=self.data[i2]))

    def computeNJTree(self):
        """comp distance matrix"""
        self.data=np.array(self._getUniqueHaplotypes().keys())
        self.nHap=self.data.shape[0]
        dMatrix=np.zeros([self.nHap,self.nHap])
        for i in np.arange(self.nHap):
            for j in np.arange(self.nHap):
                dMatrix[i,j]=self._getDistance(i,j)
        return dMatrix

#-----------------------------I/O----------------------------
    def readVCFFile(self, file, panel, selPos, seqPos, 
                    polarizeFile="/data/selectiveSweep/data/ancestral/genes_HC.txt",
                    vcfHasAA=False,excludeNonPolarizable=False):
        """read vcf data into data structure; arguments:
            file:       vcf file name
            panel:      population to use, if it is a list, or np.array
                        they will be used, otherwise determined from header
            selPos:     position of the presumably selected site
            seqPos:     start and ending position
            polarizeFile: file with info on anc/der allele
            vcfHasAA:   if Ancestral Allele is in VCF, this is used instead
            """
        polarize=False
        if not vcfHasAA:
            if polarizeFile != None:
                pf=np.loadtxt(polarizeFile,dtype="S")
                polarize=True
        if type(panel) == list or type(panel) == np.ndarray:
            individualIDs = panel
        else:
            a=np.loadtxt("interim_phase1.20101123.ALL.panel",dtype="S")
            toKeep=panel
            individualIDs=[a[i,0]  for i in range(len(a[:,1])) if a[i,1] in toKeep]

        file= open(file,"r")
        line=file.readline()
        while line[0]=="#":
            header= line
            line=file.readline()

        hd = header.split()
        data = list()
        individualsToKeep=[]
        for i,id in enumerate(individualIDs):
            for j,id2 in enumerate(hd):
                if id == id2:
                    individualstoKeep.append(j)
            

        #individualsToKeep=[i for i in range(len(hd)) if np.array(hd)[i] in individualIDs]  
        print individualsToKeep
        ls=line.split()
        ht=np.array([ ht[0:3] for ht in np.array(ls)[individualsToKeep]])
        snpdata=(np.array([i.split("|") for i in ht])).flatten()
        snpmeta=ls[:9]
        snppos=ls[1]
        snpqual=ls[6]
        snp=(snppos,snpqual,snpdata,snpmeta[3],snpmeta[4])

        data.append(snp)
        nSNPRemoved =0
        while True:
            line=file.readline()
            if not line:
                break
            ls=line.split()
            ht=np.array([ ht[0:3] for ht in np.array(ls)[individualsToKeep]])
            snpdata=(np.array([i.split("|") for i in ht])).flatten()
            snpmeta=ls[:9]
            snppos=ls[1]
            snpqual=ls[6]

            if vcfHasAA:
                moreData = ls[7].split( ";" )
                mdd = dict( [l.split("=") for l in moreData] )
                if   mdd["AA"]         == snpmeta[3]: #SNP is ancestral
                    pass
                elif mdd["AA"]         == snpmeta[4]: #SNP is derived
                    snpdata = np.array(1 - np.array(snpdata,dtype="b"),
                                       dtype="S1")
                elif mdd["AA"].upper() == snpmeta[3] and \
                        (not excludeNonPolarizable):
                    pass
                elif mdd["AA"].upper() == snpmeta[4] and \
                        (not excludeNonPolarizable):
                    snpdata = np.array(1 - np.array(snpdata,dtype="b"),
                                       dtype="S1")
                else:
                    if excludeNonPolarizable:
                        nSNPRemoved+=1
                        continue
                    #print snppos, mdd["AA"], snpmeta[3], snpmeta[4]
                    #raise ValueError("SNP neither anc nor derived")
            snp=(snppos,snpqual,snpdata,snpmeta[0],snpmeta[3],snpmeta[4])
            data.append(snp)
        print "removed %d SNP"%nSNPRemoved
        if polarize:
            chr=snpmeta[0][0]
            #refAllele = [i[3] for i in snpmeta] 
            #altAllele = [i[4] for i in snpmeta] 
            #refDict is a dict [chr, pos] => [refH, refC]
            refDict = dict()
            for i,d in enumerate(pf):
                refDict[(pf[i,0],pf[i,1])] = (pf[i,2],pf[i,3])

            toSwitch=[]
            for i,snp in enumerate(data):
                if (snp[3],snp[0]) in refDict:
                    ref = refDict[(snp[3],snp[0])]
                    #print "[%s:%s]1000k: [%s:%s] emilia: [%s:%s]" %(snp[3],snp[0],snp[4],snp[5],ref[0],ref[1])
                    if ref[0] != snp[4]:
                        pass
                        print "snp [%s:%s] is dif between 1k gen and emilia" %  (snp[3],snp[0])
                    elif ref[0] == ref[1]: #snp is polarized correctly
                        if ref[0] == snp[4]:
                            pass
                            print "snp [%s:%s] is polarized correctly" %  (snp[3],snp[0])
                        elif ref[0] == snp[4]:
                            toSwitch.append(snp[0])
                            newData = np.array(["%i"%(1-int(iii)) for iii in data[i][2]],dtype="S1")
                            t = data[i]
                            data[i] = (t[0], t[1],newData,t[3],t[4],t[5])                          
                            print "snp [%s:%s] is polarized wrong!" %  (snp[3],snp[0])
                        else: 
                            pass
                            #print "snp [%s:%s] is ODD (1)!" %  (snp[3],snp[0])
                    else: #snp is not polarized correctly
                        toSwitch.append(snp[0])
                        #print data[i]
                        print "snp [%s:%s] is polarized wrong! (2)" %  (snp[3],snp[0])
                        newData = np.array(["%i"%(1-int(iii)) for iii in data[i][2]],dtype="S1")
                        t = data[i]
                        data[i] = (t[0], t[1],newData,t[3],t[4],t[5])
                        #print data[i]
                else:
                    print "cant polarize %s:%s" % (snp[3],snp[0])

        d2=[d for d in data if (np.any(d[2]=="1") and np.any(d[2]=="0")) or int(d[0])==selPos]
        d3=d2#[d for d in d2 if d[0] == selPos or not np.all(d[2]=="1")]

            
        file.close()    


        nHap=len(d3[0][2])
        nSNP=len(d3)

        segSites=[d[0] for d in d3]
        haplotypes=list()
        for h in range(nHap): haplotypes.append([ d[2][h] for d in d3 ])
        haplotypes=np.array(haplotypes)


        self.type="vcf"
        self.nHap=nHap
        self.segSites=np.array([int(i) for i in segSites])
        self.nSegsites=nSNP
        self.data=haplotypes
        #self.individuals=np.arange(self.nHap)


        #selId=self.getClosestSNPfromPos(selPos)
        try:
            selId=self.getIDfromPos(selPos)
        except:
            raise KeyError(individualsToKeep)
        selSiteData=d3[selId][2]
        self.selId=selId
        self.selPos=selPos
        self.selSiteData=selSiteData


        self.parTheta=0
        self.parRho=0
        self.nDataSets=1
        self.length=seqPos[1]-seqPos[0]
        self.selPos-=seqPos[0]
        self.segSites=np.array([int(i)-seqPos[0] for i in self.segSites])

    def readIhsFile(self,file='ihs/sample.ihshap',snp='ihs/sample.ihsmap'):
        """reads input file as used in the ihs program by voight et al., 2006"""
        self.type="ihs"

        data=np.loadtxt(file,dtype="bool")
        snp=np.loadtxt(snp,dtype="S10")

        nHap,nSNP = data.shape
        self.snp=np.array([int(i) for i in snp[:,1]])
        snp=np.array([float(i) for i in snp[:,2]])

        
        self.data=data
        self.segSites=snp
        self.nHap=nHap
        self.nSegsites=nSNP
        #self.individuals=np.arange(self.nHap)

    def readPhaseFile(self,file='sweep-1.1/CCR5_ceu.phase',
                      snp="sweep-1.1/CCR5_ceu.snp"):
        """reads input file for sweep program in phase format"""

        self.type="phase"
        file=fileinput.input(file)


        data = np.loadtxt(file, dtype="S20")
        data = data[:,2:]
        nHap,nSNP = data.shape
        
        #convert to int
        data = [int(d) for d in data.flatten()]
        data = np.array(data)
        data=data.reshape(nHap,nSNP)

        #self.individuals=np.arange(self.nHap)

    def readNextDataSet(self):
        """reads the next data set of a multi dataset mbs output file"""
        if not hasattr(self,"file") or not hasattr(self,"nDataSets"):
            raise ValueError("no DataSet Loaded")

        if self.currDataSet >= self.nDataSets:
            raise ValueError("no more data sets left :-(")

        if not hasattr(self,"type"):
            raise ValueError("I don't know how to handle this input file type")

        if self.type=="mbs":
            self.readNextMbsDataSet()
        elif self.type=="ms":
            self.readNextMsDataSet()
        else:
            raise ValueError("Method for this input type not yet implemented")

        if hasattr(self,"pi"):
            del self.pi

        if hasattr(self,"S"):
            del self.S

        if hasattr(self,"td"):
            del self.td

        if hasattr(self,"fwh"):
            del self.fwh

        if hasattr(self,"sfs"):
            del self.sfs

        if hasattr(self,"ihs"):
            del self.ihs

        if hasattr(self,"ehh"):
            del self.ehh

        if hasattr(self,"setSFS"):
            del self.setSFS

        if hasattr(self,"freqs"):
            del self.freqs

        self.createSNPdict()
        self.currDataSet += 1

    def readNextMbsDataSet(self): 
        """get state of selected Site"
        # note that at selected Site, 1 codes for ancestral, 0 codes for derived
        # allele"""
        h1 = self.file.readline()
        self.selSiteData= [int(x) for x in
                  h1.split(':')[1].replace('d','1').replace('a','0').split()]
        self.nSegsites = int(self.file.readline().split()[1])
        "get pos of segregating sites"
        self.segSites = [int(x) for x in self.file.readline().split()[1:]]

        self.selId = sum(np.less(self.segSites, self.selPos))
        self.segSites.insert(self.selId,self.selPos)

        self.data = np.zeros((self.nHap,self.nSegsites),dtype="bool")
        #self.data[:,0]=[int(i) for i in self.segSites]
        if len(np.unique(self.segSites)) != len(self.segSites):
            #print "correcting",
            for i,s in enumerate(self.segSites):
                if sum(s==np.array(self.segSites))>1:
                    #print self.segSites[i]
                    self.segSites[i+1]+=1
        c=0
        while 1:
            buf=self.file.readline()
            if buf=="\n" or not buf:
                break
            buf = buf.strip()
            s = [int(buf[i]) for i in range(len(buf))]
            s.insert(self.selId,self.selSiteData[c])
            self.data[c,:]=s
            c=c+1

    def readMbsFile(self,file):
        """reads output from MBS"""
        self.file = open(file)
        self.type = "mbs"

        #first, parse the command as argument
        h1=self.file.readline()
        h1=h1.split()
        tpos,fpos,spos,rpos =(-1,-1,-1,-1)

        for index, value in enumerate(h1):
            if value=='-t': tpos=index
            if value=='-f': fpos=index
            if value=='-s': spos=index
            if value=='-r': rpos=index

        self.parTheta = float(h1[tpos+1])
        self.parRho = float(h1[rpos+1])
        self.nDataSets = int(float(h1[fpos+1]))
        
        self.currDataSet = 0 
        self.selPos = int(float(h1[spos+2]))
        self.length = int(float(h1[spos+1]))
        self.nHap= int(float(h1[2]))
        #self.individuals=np.arange(self.nHap)
	self.trajName=h1[spos+3]
	if self.trajName[0:3]=="chr": 
	    self.chromosome=self.trajName[3:]
	    #print "inferred chr:", self.chromosome


        ln=self.file.readline()
        if ln!="\n":
            print ln
            #assume it is individual ids
            self.individualIDs=ln.split(",")

        #read first DataSet
        self.readNextDataSet()
        
    def readMsFileCombine(self,file):
        self.file   =   open(file)
        self.type   =   "msc"
        self.length =   1
        h1          =   self.file.readline()
        h1          =   h1.split()

        self.nHap   =   int(float(h1[1]))
        self.nDataSets= 1
        self.nSegsites= int(float(h1[2]))
        self.data=np.zeros((self.nHap,self.nSegsites),dtype="b")
        self.segSites   = np.arange(self.nSegsites,dtype="f")/self.nSegsites

        i = 0 
        while i<self.nSegsites:
            buf=self.file.readline()
            if buf =="":
                break
            if buf[0:8]=="segsites":
                self.file.readline()
                c=0
                while 1:
                    buf=self.file.readline()
                    if buf=="\n" or not buf:
                        break
                    buf = buf.strip()
                    s = int(buf)
                    self.data[c,i]=s
                    c+=1
                i+=1



    def readMsFile(self,file):
        self.file = open(file)
        self.type = "ms"
        self.length=1

        #parse the command line:
        h1=self.file.readline()
        h1=h1.split()

        self.nHap = int(float(h1[1]))
        self.nDataSets = int(float(h1[2]))
        self.currDataSet = 0
        #self.individuals=np.arange(self.nHap)

        self.readNextDataSet()

    def readNextMsDataSet(self):
        if hasattr(self,"data"):
            del self.data
        if hasattr(self,"segSites"):
            del self.segSites
        if hasattr(self,"nSegSites"):
            del self.nSegsites
        if hasattr(self,"SNPdict"):
            del self.SNPdict
        if hasattr(self,"selId"):
            del self.selId
        if hasattr(self,"selPos"):
            del self.selPos
        if hasattr(self,"selSiteData"):
            del self.selSiteData
        while 1:
            buf=self.file.readline()
            if buf[0:8]=="segsites":
                buf=buf.split()
                self.nSegsites=int(float(buf[1]))
                print self.nSegsites
                break
            if buf == '':
                raise ValueError("could not read input file: number of\
                                     segsites not found")

        #read the positions
        self.segSites = [float(x) for x in self.file.readline().split()[1:]]
        #check if all segregating sites are unique, otherwise it will not work:
        if len(np.unique(self.segSites)) != len(self.segSites):
               raise ValueError("some segregating sites have the same index//ms")

        #prepare Data
        self.data = np.zeros((self.nHap,self.nSegsites),dtype="int")

        c=0
        while 1:
            buf=self.file.readline()
            if buf=="\n" or not buf:
                break

            buf = buf.strip()
            s = [int(buf[i]) for i in range(len(buf))]
            print len(s)
            self.data[c,:]=s
            c=c+1

        self.setSelectedSite(pos=.5)

    def readNeurospora(self,file="Ellison_2011_SNPdata.txt"):
        self.type="neurospora"
        data=[]
        pos=[]
        chr=[]
        def toInt(c):
            return int(c)
        
        data = np.loadtxt(file,dtype="S20",comments="#", converters = {0:toInt,
                                                                       1:toInt})

        self.chrSNP=data[:,0]
        self.segSites=data[:,1]
        data=data[:,2:]
        self.nSegsites,self.nHap=data.shape
        major=np.zeros(self.nSegsites)
        minor=np.zeros(self.nSegsites)
        self.data=np.zeros(data.shape,dtype="i1")
        for i in range(self.nSegsites):
            #self.data[i],b,(major[i],minor[i])    = self.makeNumeric(data[i])
            self.data[i],b,_    = self.makeNumeric(data[i])
            if not b: print i

        self.data=np.transpose(self.data)



    def readPhasedFile(self,file='phased/AllDiAllelicSitesEpas1.tibet.like.phased',selPos=0,filter=False):
        """reads a file from the (phased) epas 1 gene, by beagle program"""
        self.type="beagle"

        #get chromosome information
        f=open(file,"r")
        s=f.readline()
        s=f.readline()
        ss=s.split()
        self.chromosome=ss[1].split("_")[0][3:]
        f.close()
        #print "chromosome inferred: ", self.chromosome
        file=fileinput.input(file)

        header = file.readline()
        header = header.split()

        def id2pos(id):
            return int(id.split('_')[1])

        #gets the individuals
        for i in range(len(header[3::2])):
            header[2+i*2] += '_a'
            header[3+i*2] += '_b'

        header[1]='pos'

        self.header = header
        self.individualIDs = self.header[2:]

        #annotate samples to pops
        self.han=np.where([s[:2]=="NA" for s in self.individualIDs])[0]
        self.tib =np.where([s[:2]=="DR" or s[:2]=="NQ" for s in
                                  self.individualIDs])[0]
        self.pops=np.zeros(len(self.individualIDs))
        self.pops[self.tib]=1


        type=['S1','I16']; type.extend(['S1' for x in range(len(header)) if x > 1])

        dt={'names':header,'formats':type}
        data = np.loadtxt(file, dtype="S10" ,
                          comments='I', converters = {1:id2pos} )




        #shift to start from 0/x offset
        pos = data[:,1]
        pos = np.asarray([ int(s) for s in pos ])
        #self.minPos = 46358067
        self.minPos = 0
        pos = pos - self.minPos
        self.length = max(pos)

        data= data[:,2:]
        nSNP,nHap = data.shape
        #self.individuals=np.arange(nHap)
        data_numeric = np.zeros((nSNP,nHap),dtype='uint8')



        badData = np.zeros(nSNP,dtype="bool")
        self.major = np.empty(nSNP,dtype="S1")
        self.minor = np.empty(nSNP,dtype="S1")
        for i in range(nSNP):
            data_numeric[i,:],badData[i],(self.major[i],self.minor[i]) = makeNumeric(data[i,:])
            if self.verbose: print(i,badData[i])

        file.close()

        #remove all monomorphic SNP
        if filter:
            data_numeric = data_numeric[badData]
            pos = pos[badData]

        self.data = data_numeric.transpose()
        self.nHap,self.nSegsites = self.data.shape

        self.segSites = pos
        self.createSNPdict()
         #use arbitrary SNP pos for now
        if selPos==0:
            self.setSelectedSite(pos=min(self.segSites))
        else:
            self.setSelectedSite(pos=selPos)

    def makeNumeric(self,row):
        nucleotides=np.array(('A','C','G','T'))
        b=1
        alleleFreq=np.array(((sum(row=='A'),sum(row=='C'),sum(row=='G'),sum(row=='T'))))
        if(sum(alleleFreq != 0) != 2):
            if self.verbose: print 'WAAAAAAAAAAAAAAAH: Bad genotype',alleleFreq
            b=0
        k=np.argsort(alleleFreq)[::-1]
        major=nucleotides[k[0]]
        minor=nucleotides[k[1]]
        if alleleFreq[k[1]]==0:
            minor="N"
        nuc=nucleotides[k]
        return (np.where(row==major,0,1),b,(major,minor))

    def writeToHaploViewFile(self,file,data=None):
        if data == None:
            data = self.data

        nHap,nPos = data.shape
        families = ['f'+str(i/2) for i in range(nHap)]
        ids = ['i'+str(i/2) for i in range(nHap)]

        output = np.empty((nHap,nPos+2),dtype="S10")
        output[:,0] = np.asarray(families)
        output[:,1] = np.asarray(ids )
        output[:,2:] = data+1


        np.savetxt(file,output,fmt="%s")

    def writeToMsFile(self,file=None,data=None,segSites=None,par=None,selSiteInData=True):
        """if selSiteInData, the selected Site is writen into the data portion
        as well"""
        if data == None:
            data = self.data
        nHap,nPos = data.shape
        if par == None:
            dummy="dummy"
            if hasattr(self,"chr"): dummy="chr"+self.chr
            par=(nHap,self.parTheta,self.parRho,self.nDataSets,self.length,self.selPos,dummy)
        if segSites == None: 
            segSites = self.segSites

        nSNP=len(segSites)


        selSite="".join((str(i)+" " for i in self.selSiteData))
        selSite=selSite.replace("1","d").replace("0","a")

        if file == None:
            #first line fake ms command
            print "./ms %i 1 -t %f -r %f -f 1 %i -s %i %i %s" % par
            print "110 110 110"
            print "\n"
            print "//"
            print "segsites: %i" % nPos
            s="positions: ";
            if segSites != None:
                for ss in segSites:
                    if ss!=self.selPos or selSiteInData:
                        s+=str(ss)+" "
                
            print s 
            for i in range(nHap):
                if selSiteInData:
                    s = "".join([str(int(j)) for j in data[i,:]])
                else:
                    s = "".join([str(int(j)) for j in data[i,:self.selId]])
                    s += "".join([str(int(j)) for j in data[i,self.selId+1:]])
                print s
        else:
            out=open(file,'w')
            #first line fake ms command
            out.write("./ms %i 1 -t %f -r %f -f 1 %i -s %i %i %s\n" %
                      par)
            out.write("110 110 110\n")
            out.write("\n")
            out.write("//\n")
            out.write("segsites: %i\n" % nPos)
            s="positions: ";
            if segSites != None:
                for ss in segSites:
                    if ss!=self.selPos or selSiteInData:
                        s+=str(ss)+" "
                
            out.write(s+"\n") 
            for i in range(nHap):
                if selSiteInData:
                    s = "".join([str(int(j)) for j in data[i,:]])
                else:
                    s = "".join([str(int(j)) for j in data[i,:self.selId]])
                    s += "".join([str(int(j)) for j in data[i,self.selId+1:]])
                out.write(s+"\n") 
    def writeToMbsFile(self,file=None,data=None,segSites=None,par=None,selSiteInData=False):
        """if selSiteInData, the selected Site is writen into the data portion
        as well"""
        if data == None:
            data = self.data
        nHap,nPos = data.shape
        if par == None:
            dummy="dummy"
            if hasattr(self,"chr"): dummy="chr"+self.chr
            par=(nHap,self.parTheta,self.parRho,self.nDataSets,self.length,self.selPos,dummy)
        if segSites == None: 
            segSites = self.segSites

        nSNP=len(segSites)


        selSite="".join((str(i)+" " for i in self.selSiteData))
        selSite=selSite.replace("1","d").replace("0","a")

        if file == None:
            #first line fake ms command
            print "command:\t./mbs %i -t %f -r %f -f 1 %i -s %i %i %s" % par
            if hasattr(self,"individualIDs"):
                print(",".join(self.individualIDs))
            else:
                print("")
            print "//0-1 allels: "+selSite
            print "segsites: %i" % nPos
            s="positions: ";
            if segSites != None:
                for ss in segSites:
                    if ss!=self.selPos or selSiteInData:
                        s+=str(ss)+" "
                
            print s 
            for i in range(nHap):
                if selSiteInData:
                    s = "".join([str(int(j)) for j in data[i,:]])
                else:
                    s = "".join([str(int(j)) for j in data[i,:self.selId]])
                    s += "".join([str(int(j)) for j in data[i,self.selId+1:]])
                print s
        else:
            out=open(file,'w')
            #first line fake ms command
            out.write("command:\t./mbs %i -t %f -r %f -f 1 %i -s %i %i %s\n" %
                      par)
            if hasattr(self,"individualIDs"):
                out.write(",".join(self.individualIDs))
            else:
                out.write("\n")
            out.write("//0-1 allels: "+selSite+"\n")
            out.write("segsites: %i\n" % nPos)
            s="positions: ";
            if segSites != None:
                for ss in segSites:
                    if ss!=self.selPos or selSiteInData:
                        s+=str(ss)+" "
                
            out.write(s+"\n") 
            for i in range(nHap):
                if selSiteInData:
                    s = "".join([str(int(j)) for j in data[i,:]])
                else:
                    s = "".join([str(int(j)) for j in data[i,:self.selId]])
                    s += "".join([str(int(j)) for j in data[i,self.selId+1:]])
                out.write(s+"\n") 

    def writeSweepFile(self,file='tibet',chr=2):
        #write .snp file first, which has 3cols: snpid, chr, pos
        output = np.empty((self.nSegsites+1,3),dtype="S20")
        snpid=["snp_"+str(i) for i in range(self.nSegsites)]
        output[0:] = np.asarray(("snpid","chr","HG16"))
        output[1:,0] = np.asarray(snpid)
        output[1:,1] = np.repeat(chr,self.nSegsites)
        output[1:,2] = np.asarray(self.segSites,dtype="S20")
        np.savetxt(file+".snp",output,fmt="%s", delimiter="\t")
        
        #write data file
        output = np.empty((self.nHap,self.nSegsites+2),dtype="S20")
        indid=["ind_"+str(i) for i in range(self.nHap/2)]
        indid=np.repeat(indid,2)
        output[:,0]=indid
        output[:,1]=np.tile(("T","U"),self.nHap/2)
        output[:,2:]=self.data+1
        np.savetxt(file+".phase",output,fmt="%s", delimiter="\t")

    def writeSweepFinderFile(self,file="sweepfinder_input.txt"):
        output=np.empty((self.nSegsites+1,4),dtype="S10")
        output[0,:]=("position","x","n","folded")
        output[1:,0]=self.segSites
        self.createSFS()
        output[1:,1]=[int(i) for i in self.freqs]
        output[1:,2]=np.repeat(self.nHap,self.nSegsites)
        output[1:,3]=np.repeat(1,self.nSegsites)
        np.savetxt(file,output,fmt="%s",delimiter="\t")

    def writeArpFile(self,file):
        out=open(file,'w')
        out.write( """
        [Profile]
        Title="auto"
        NbSamples=2
        GenotypicData=0
        GameticPhase=1
        DataType=STANDARD
        LocusSeparator=WHITESPACE
        MissingData='?'

        """)
        out.write( """
        [Data]
        [[Samples]]
        SampleName="Tibet"
        """)
        out.write( "SampleSize=%i\n" % (len(self.tib),))
        out.write( "SampleData={\n")

        for i in self.tib:
            out.write( "%i %i %s" %(i, 1, " ".join([str(int(s)) for s in
                                                     self.data[i,:]])))
            out.write("\n")
        out.write( "\n}\n")

        out.write( "SampleName=\"Han\"\n")
        out.write( "SampleSize=%i\n" % (len(self.han),))
        out.write( "SampleData={\n")

        for i in self.han:
            out.write( "%i %i %s" %(i, 1, " ".join([str(int(s)) for s in
                                                     self.data[i,:]])))
            out.write("\n")
        out.write( "\n}\n")

    def writeNetworkFile(self):
        print "  ;1.0"
        print " ;".join([str(int(s)) for s in self.segSites])
        print ";".join(["10" for s in self.segSites])

        for i in self.tib:
            print ">tib_%i ;10;0;;;;;;" % i
            print "".join([str(int(s)) for s in self.data[i,:]])
        for i in self.han:
            print ">han_%i ;10;0;;;;;;" % i
            print "".join([str(int(s)) for s in self.data[i,:]])

    def readRefSeqEmilia(self,file):
        rawData=np.loadtxt(file,dtype="S")
        pos         =   np.array([int(i) for i in rawData[:,1]])
        refRaw      =   rawData[:,2]
        for i,d in enumerate(refRaw):
            if d=="NA":
                refRaw[i]=False
        a=dict([(j,refRaw[i])for i,j in enumerate(pos)])
        for i in a:
            if a[i]=="False":
                a[i]=False
        self.ancestralAllel=a


    def readRefSeq(self,file="EPAS1_BGIsnp_orthNuc"):
        rawData=np.loadtxt(file,dtype="S")
        pos         =   np.array([int(i) for i in rawData[1:,1]])
        refHum      =   rawData[1:,2]
        refChimp    =   rawData[1:,3]
        refMaq      =   rawData[1:,4]
        
        polarizable =   refChimp==refMaq
         
        ancestralAllel=dict()
        for i in np.arange(len(pos)):
            if polarizable[i] and refChimp[i] != '-':
                ancestralAllel[pos[i]]=refChimp[i]
            else:
                ancestralAllel[pos[i]]=False
        self.ancestralAllel=ancestralAllel

    def polarizeData(self,keepNA=True,verbose=True):
        if not hasattr(self,"ancestralAllel"):
            raise ValueError("no Data available. Use readRefSeq() to load a file")

        aa=self.ancestralAllel

        toRemove=[]
        self.polarizable=np.zeros(self.nSegsites,dtype="Bool")
        for i in np.arange(len(self.segSites)):
            snp=self.segSites[i]
            if not snp in aa:
                if verbose: print "SNP %i(%i) has no polarization info:\
                   [m:%s]"\
                %(i,self.segSites[i],self.major[i])
                if not keepNA:
                    toRemove.append(i)
                continue
            if aa[snp] == False:
                if verbose: print "SNP %i(%i) is missing data:\
                   [m:%s],[a:%s]"\
                %(i,self.segSites[i],self.major[i],aa[snp])
                if sum(self.data[:,i])>self.nHap/2:
                    self.data[:,i] = 1-self.data[:,i]
                    print sum(self.data[:,i])
                if not keepNA:
                    toRemove.append(i)
                continue
            if self.major[i] == aa[snp]:
                if verbose: print "SNP %i(%i) is polarized as expected:\
                   [m:%s],[a:%s]"\
                %(i,self.segSites[i],self.major[i],aa[snp])
                self.polarizable[i] = True
            elif self.minor[i] == aa[snp]:
                if verbose: print "SNP %i(%i) is polarized differently:\
                   [m:%s],[a:%s]"\
                %(i,self.segSites[i],self.major[i],aa[snp])
                self.polarizable[i] = True
                
                #adjust data
                self.data[:,i] = 1-self.data[:,i]
            else:
                if verbose: print "SNP %i(%i) is different:\
                   [maj:%s],[min:%s],[m[a:%s]"\
                %(i,self.segSites[i],self.major[i],self.minor[i],aa[snp])
                if not keepNA:
                    toRemove.append(i)
        for i in toRemove[::-1]:
            self.removeSNP(i)

    def polarizeFromDBSNP(self,f, offset=0):
        """f should be a list of positions of SNP to switch, nothing else
           offset is the starting point if the SNP are in genomic and not local
           coordinates"""
        positions=np.loadtxt(f)
        for i in positions:
            if not i-offset in self.segSites:
                continue
            id=self.getIDfromPos(i-offset)
            self.data[:,id] = 1 - self.data[:,id]
            print "SNP " ,id, "/",i,"/",i-offset, " adjusted"

    def simplifyIndividualIDs(self):
        """function to simplify individualIDs to remove the annoying beagle
        attachments"""
        self.individualIDs=[s.split(".")[0] for s in self.individualIDs]
        self.individualIDs=[s.replace("-",".") for s in self.individualIDs]
        self.individualDict=dict([(s,((i-1),i)) for i,s in enumerate(self.individualIDs)])

    def getLDStats(self,ids):
        data=np.empty((self.nHap,2),dtype="bool")

        data[:,0]=self.getDataFromId(ids[0])
        data[:,1]=self.getDataFromId(ids[1])

        n=float(self.nHap)
        p1=sum(data[:,0])/n
        p2=1-p1
        q1=sum(data[:,1])/n
        q2=1-q1
        p11=sum(np.logical_and(data[:,0]==1,data[:,1]==1))/n
        p10=sum(np.logical_and(data[:,0]==1,data[:,1]==0))/n
        p01=sum(np.logical_and(data[:,0]==0,data[:,1]==1))/n
        p00=sum(np.logical_and(data[:,0]==0,data[:,1]==0))/n

        D=p11-p1*q1
        if D<0:
            Dprime=D/np.min(p1*q1,p2*q2)
        elif D>0:
            Dprime=D/np.min(p1*q2,p2*q1)
        else:
            Dprime=0.0

        r=D/np.sqrt(p1*p2*q1*q2)
        return D,Dprime,r

    def readArpFile(self,file="sGeneticsOutput/_spc_samp_1.arp"):
        f       =   open(file)
        line    =   f.readline()
        data    =   []
        sampleId=[]
        while (line !=""):
            if line[0]=="\n" or line[0]=="#":
                line=f.readline()
                continue
            if "SampleName" in line:
                sampleId.append(line.split("\"")[1])
            if "SampleData" in line:
                data.append(f.readline().split()[2:])
                data.append(f.readline().split())
            line=f.readline()
        f.close()

        self.data=np.array(data,dtype="i")
        self.nHap,self.nSegsites=self.data.shape


    def addRecombinationMap(self,startPos=39586954,
                             file="/data/selectiveSweep/IFNL/DATA/genetic_map_GRCh37_chr19.txt"):
        """
        calculates the position of each SNP on a recombination map. This is done
        very lazily in an O(n*m) algorithm, where n is the number of map entries
        and m is the number of SNP. Should be improved when this becomes time
        critical
        """
        geneticMap=[]
        f = open(file)
        for line in f:
            line = line.split()
            curLine=[line[1],line[4]]
            geneticMap.append(curLine)

        geneticMap=np.array(geneticMap,dtype="f4")
        
        segSites_Rec=[]
        for c,curPos_Physical in enumerate(self.segSites):
            print c,"/",self.nSegsites
            i=0
            #lazily ignoring the case where the data extends over the chromosome
            while curPos_Physical+startPos > geneticMap[i][0]:
                i+=1
            frac = (curPos_Physical+startPos - geneticMap[i-1][0]) / \
                   (geneticMap[i][0] - geneticMap[i-1][0])
            curPos_Rec = geneticMap[i-1] + frac * (geneticMap[i] - geneticMap[i-1])
            segSites_Rec.append(curPos_Rec[1])
            
            

