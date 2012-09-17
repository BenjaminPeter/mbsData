#!/usr/bin/env python

from mbsData import mbsData
from mbsData2pop import mbsData2P

class mbsData2PSSS(mbsData2P):
    """this should allow usage of the standard mbsData function using the output from the SSSimulator. Currently, data comes in the following form:
        1. 2D SFS, ending *.sfs
            this files have a 2dSFS, the first line is the total tree length, the second through nth lines are the actual SFS
        2. Frequency tables, ending *.ft
            each row of a FT is a tree segment, giving the lenght of the tree segment followed by the (absolute) allele frequency in each population. Sample Sizes have to be given externally...
        3. output of statistics, ending *.stat
            gives popid1, popid2 FST sdFST, psi sdpsi, deltaH, sddeltahfor all pairwise comparisons

10.9.2012: currently solely working on FT
    """
    def readFT(self,file="out__1.ft"):
        """reads file format where the first column is the branch length, following columns are the absolute frequency in all populations. Currently sample sizes in all pops are assumed to be 100"""



