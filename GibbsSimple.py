#By Ben Battenburg
#Computational Genomics - SimpleGibbs Implementation
#This program detects an over-present motif using a W (width) parameter of a predicted motif in a seqdata.txt file and a background pseudocount parameter. Sequence data is formatted as one sequence per line.

#Notes: The longer the motif width, the greater phase shift problem (stuck in local optima). 
#A smaller background parameter causes less fluctuations in random indexes, i.e., it converges faster but I assume is less specific in detection. 
#Larger background parameters seem to decrease the likelihood for convergence. Background parameter seems to limit the phase shift problem when the right value is applied.
#The larger the sequences, the smaller the background parameter in order for convergence.

import random
import math
import bisect
from operator import mul

#Parameters main
wMTheta = 10 #Motif width
bgTheta = 1.3 #Background psuedocount
##rtTheta = 400 #Number of iterations in a run. - Replaced by autoconverge.
cVTheta = 50 #Number of most recent past iterations with same output to tell convergence.

#Parameters for Phase Shift correction.
mSTheta = 10 #Number of motifs sampled on each side of current motif.
pSTheta = 0 #Starting iteration to begin phase shift.
pFTheta = 5000 #Final iteration to do phase shift
pITheta = 20 #Iteration Hz for shift between pSTheta and pFTheta. *Must be > 0
sMTheta = 2 #Positive score for matched residue motif position.
sNTheta = -1 #Negative score for mismatched residue motif position.

#Debug Parameters
genMotif = ['GATTACATAC','GATTACATAC','GATTACATAC','GATTACATAC']
testList = []
testInt = 0
#Generated Motif indices for each sequence.
#287  
#235  
#370 
#53    

class tableCalc():
    def __init__(self, seqDataExcl,seqData,seqPosR,seqRandIndex,seqCount, randIndexExcl,wMTheta,bgTheta):
        #Input class variables
        self._seqDataExcl = seqDataExcl #Excludes current sequence
        self._randIndexExcl = randIndexExcl #Excludes current motif position index
        self._seqDataCurrent = seqData #Current sequence in string format
        self._seqDataCurrentP = seqPosR #Random initialized start position for motif
        self._seqCount = seqCount #Current sequence index from seqData
        self._seqRandIDCur = seqRandIndex #Current sequence motif index
        self._wMTheta = wMTheta #Motif width parameter
        self._bgTheta = bgTheta #Background parameter
        
        #Class specific variables. - Yikes...
        self._res = ['A','C','G','T']
        self._bkgrtotA = 0
        self._bkgrtotC = 0
        self._bkgrtotG = 0
        self._bkgrtotT = 0
        self._bkgrA = 0
        self._bkgrC = 0
        self._bkgrG = 0
        self._bkgrT = 0    
        self._bkgrListACGT = []
        self._count = 0
        self._table = []
        self._seqPosR = []
        self._tableRow = []
        self._rowA = []
        self._rowC = []
        self._rowG = []
        self._rowT = []
        self._tableF = []
        self._predUpdTemp = []
        self._predUpd = []
        self._predUpdBg = []
        self._SamPosList = []
        self._seqWeights = []
        self._totals = []
      
    #Create background profile utilizing bgTheta.   
    def bkground(self):
    #Counts the total number of each residue in sequence and subtracts from it the occurence of each residue in the purposed motif.
        #Slice purposed motifs
        for n in self._seqDataExcl:
            self._seqPosR.append(n[self._randIndexExcl[self._count]:self._randIndexExcl[self._count]+self._wMTheta])
            self._count += 1         
        #Count residues in motif.
        for n in self._seqPosR:
            self._bkgrA += n.count('A')
            self._bkgrC += n.count('C')
            self._bkgrG += n.count('G')
            self._bkgrT += n.count('T')        
        #Count total residues.
        for n in self._seqDataExcl:
            self._bkgrtotA += n.count('A')
            self._bkgrtotC += n.count('C')
            self._bkgrtotG += n.count('G')
            self._bkgrtotT += n.count('T')
            
        #Calculate background list.
        self._bkgrListACGT.append(self._bkgrtotA-self._bkgrA)
        self._bkgrListACGT.append(self._bkgrtotC-self._bkgrC)
        self._bkgrListACGT.append(self._bkgrtotG-self._bkgrG)
        self._bkgrListACGT.append(self._bkgrtotT-self._bkgrT)
    
    #Make tables for residue frequency counts, predictive updates and sampling    
    def probTable(self):
               #Probability table
        #Creates a list of individual motif residues and appends each motif residue list to another list -> _table.
        for n in self._seqPosR:
            self._listtemp = list(n) 
            self._table.append(self._listtemp)
        
        #Reformat motifs to column format for residue frequency count.
        self._cols = map(list, zip(*self._table))    
        #Counts residues by parsing A,C,G,T by column and appends each count to a single list (e.g., 1,0,1,1 for the first four appended list elements means there is 1 A, 0 C, 1 G, and 1 T in _cols[0]).
        for m in range(len(self._cols)):
            for n in self._res:       
                self._tableRow.append(self._cols[m].count(n))
        #Separates residues from _tableRow by intervals correlating to specific appended residues (In this case 4 for A,C,G,T) and creates a list for each residue count.       
        self._rowA = self._tableRow[0:len(self._tableRow):4]
        self._rowC = self._tableRow[1:len(self._tableRow):4]
        self._rowG = self._tableRow[2:len(self._tableRow):4]
        self._rowT = self._tableRow[3:len(self._tableRow):4]
        #Appends each residue row to a Final table.
        self._tableF.append(self._rowA)
        self._tableF.append(self._rowC)
        self._tableF.append(self._rowG)
        self._tableF.append(self._rowT)
        
        #Calculate residue motif position profile for predictive update by looping through _tableF and using background parameters.        
        for j in range(len(self._tableF)):
            for i in range(len(self._tableF[j])):
                self._predUpdTemp.append((self._tableF[j][i] + (self._bgTheta/4.0)) / ((len(self._seqDataExcl)) + self._bgTheta))
            self._predUpd.append(self._predUpdTemp)
            self._predUpdTemp = []
        
        #Calculate background profile.
        self._count = 0 #Count for background total.
        
             #Calculate background total for pseudocount calculation.
        for i in self._bkgrListACGT:
            self._count += i 
             #Pseudo count calculation/normalization.
        for j in range(len(self._bkgrListACGT)):
            self._predUpdBg.append((self._bkgrListACGT[j] + self._bgTheta/4.0) / (self._count + self._bgTheta))
        
    #Sample motif positions, determine weights, and randomly select next position.  
    def samplePos(self):
        self._sListtemp = []
        self._posList = []
        #Construct list of possible motif positions in sequence of interest.
        for j in range((len(self._seqDataCurrent) - self._wMTheta) + 1):
            self._posList.append(self._seqDataCurrent[j:j+self._wMTheta])
            continue
            
         #Parse through each motif and calculate position weights.              
        for s in self._posList:   
            for i in range(self._wMTheta):       
                if s[i] == 'A':
                    self._sListtemp.append(self._predUpd[0][i] / self._predUpdBg[0])
                elif s[i] == 'C':
                    self._sListtemp.append(self._predUpd[1][i] / self._predUpdBg[1])              
                elif s[i] == 'G':
                    self._sListtemp.append(self._predUpd[2][i] / self._predUpdBg[2])                           
                elif s[i] == 'T':
                    self._sListtemp.append(self._predUpd[3][i] / self._predUpdBg[3])              
                #Construct Position weights.
            self._SamPosList.append(reduce(mul,self._sListtemp,1))
                #Reset sListtemp
            self._sListtemp = []
            
        #Normalize motif position weights and randomly select new position using binary search.
             #Normalization
        self._seqWeights = [x/sum(self._SamPosList) for x in self._SamPosList]
        runningTotal = 0
             #Binary search selection.    
        for t in self._seqWeights:
            runningTotal += t
            self._totals.append(runningTotal)
        binRan = random.random()
        return bisect.bisect_right(self._totals,binRan)
    
    #Motif phase shift to avoid local optima.
    def phaseShift(self, randIndex):
        self._randIndex = randIndex
        self._randIndexShiftL = []
        self._randIndexShiftR = [] 
        self._scoreSheet = []
        self._shiftIndices = []
        self._count = 0        
        
        #Build alternate shifted motifs.
        for m in range(1,mSTheta+1):
            self._randIndexShiftL.append([x-m for x in self._randIndex])
            self._randIndexShiftR.append([x+m for x in self._randIndex])
            
        #Test max and min values before grading. - Avoids index out of range.
        for s in [self._randIndexShiftL, self._randIndexShiftR]:
            for t in s:
                for r in t:
                    if r < 0 or r > len(self._seqDataCurrent) - self._wMTheta:
                        del s[self._count:]
                        break
                    else:
                        continue
                self._count += 1
            self._count = 0
            
        for s in [self._randIndexShiftL, self._randIndexShiftR]:
            if not s:
                continue
            else:
                for t in s:
                    self._scoreSheet.append(tblCalc.phaseScore(t, seqData))
        #Check if scoreSheet has values to compare.
        if not self._scoreSheet:
            return self._randIndex
        else:
            #Find best motif from scoreSheet and compare it to the current position score.
            self._shiftIndices = self._randIndexShiftL + self._randIndexShiftR
            if max(self._scoreSheet) >= tblCalc.phaseScore(self._randIndex, seqData):
                return self._shiftIndices[self._scoreSheet.index(max(self._scoreSheet))]
            else:
                return self._randIndex
        
    #Grade phase shifted motifs and return a total score.    
    def phaseScore(self, randIndex, seqData):
        self._Indices = randIndex
        self._seqData = seqData
        self._mPosList = []
        self._curList = []
        self._excList = []
        self._count = 0
        self._mScore = 0  
        self._n = 0  
        
        for m in self._seqData:
            self._mPosList.append(m[self._Indices[self._n]:self._Indices[self._n]+wMTheta])
            self._n += 1
        
        for l in range(len(self._mPosList)):
            self._curList = list(self._mPosList[l])    
            self._excList = [list(x) for x in self._mPosList]
            del self._excList[l]
            
            for a in self._curList:
                for s in range(len(self._excList)):
                    if a == self._excList[s][self._count]:
                        self._mScore += sMTheta        
                    else:
                        self._mScore -= sNTheta
                        
                self._count += 1
            self._count = 0
            
        #Makes sure score is not zero (Possible error when comparing motifs and null lists)
        if self._mScore == 0:
            self._mScore +=1
            return self._mScore
        else:
            return self._mScore     
        
    #Used in debugging parameters
    def scorePrint(self, finalM):
        self._finalM = finalM
        self._curList = []    
        self._excList = []
        self._count = 0 
        self._mScore = 0
        
        for l in range(len(self._finalM)):
            self._curList = list(self._finalM[l])    
            self._excList = [list(x) for x in self._finalM]
            del self._excList[l]
                    
            for a in self._curList:
                for s in range(len(self._excList)):
                    if a == self._excList[s][self._count]:
                        self._mScore += sMTheta                              
                    else:
                        self._mScore -= sNTheta
                                
                self._count += 1
            self._count = 0                     
        return self._mScore
    
    #Automate convergence check.
    def convgCheck(self, randIndex,conLog):
        self._conLog = conLog
        self._randIndex = randIndex
        #Make list to log and count past cVTheta motif indices.
        self._conLog.append(self._randIndex[:])
        #Count        
        print self._conLog.count(self._randIndex)
        if self._conLog[-cVTheta:].count(self._randIndex) >= cVTheta:
            return True
        else:
            return False
#
                                 ##Main
#
#Open sequence data from seqData.txt.
f = open('seqData.txt', 'r')
#Initialize variables.
seqData = []
seqDataExcl = []
seqRandIndex = 0
ConVg = False
Motifs = []
conLog = []

#Get line sequence data from seqData.txt file(i.e, 1 line = 1 sequence)
for line in f:   
    seqData.append(line.rstrip())

                  #Initialize sample positions.
#Select random position for each possible motif having width wMTheta within each sequence.
randI = 0 #random integer for choosing beginning motif index.
randIndex = [] #List of beginning motif indices in sequence.
seqPosR = [] #Compiled list of random motif positions.

for s in seqData:
    #"len(s) - wMTheta": Number of possible starting motif indices within sequence.
    randI = random.randint(0, (len(s) - wMTheta)) 
    randIndex.append(randI)    
    seqPosR.append(s[randI:randI+wMTheta])
    
#
          #Begin loop and construct class for sampling.           
#
#Variable initializtion
seqCount = 0 #Specific sequence currently being iterated.
term = 0 #Debug iteration count.

while ConVg == False:

    #Save current sequence motif start position for update.
    seqRandIndex = randIndex[seqCount]
    #Exclude current motif start position from list to pass into class.
    randIndexExcl = list(randIndex)
    del randIndexExcl[seqCount]    
    #Exclude current sequence from list to pass into class.
    seqDataExcl = list(seqData) 
    del seqDataExcl[seqCount]
    
    #Construct class to calculate probabilities and to sample/update new motif position for the specific sequence iteration.
    tblCalc = tableCalc(seqDataExcl,seqData[seqCount],seqPosR[seqCount],seqRandIndex,seqCount,randIndexExcl,wMTheta,bgTheta)
    
    tblCalc.bkground() #Calculates background profile.
    tblCalc.probTable() #Constructs tables for sampling.
    
    #Samples each possible position for current sequence motif, creates weights, and randomly returns new motif start position.
    randIndex[seqCount] = tblCalc.samplePos() 
        
    #Calling phaseShift to avoid local optima.
    a = range(pSTheta,pFTheta+pITheta,pITheta)
    for i in a:
        if term == i:
            randIndex = tblCalc.phaseShift(randIndex)
        else:
            continue    

    print 'i:' + str(term)          
    print randIndex 

    #Check last cVTheta iterations for convergence.
    ConVg = tblCalc.convgCheck(randIndex,conLog)     
    
    #Iterative conditionals for going through sequential sequences.
    if seqCount == len(seqData) - 1:
        seqCount = 0      
    else:
        seqCount += 1
    
    if ConVg == True: 
        for s in range(len(randIndex)):
            Motifs.append(seqData[s][randIndex[s]:randIndex[s]+wMTheta])
        #Print found motifs.
        for p in range(len(randIndex)):
            print 'Sequence:' + str(p+1) + ', Residue:' + str(randIndex[p]) + ', Motif:' + str(Motifs[p])
        #Debug code    
        print '\nOutput score: ' + str(tblCalc.scorePrint(Motifs)) + '.'
        print 'Generated motif score: ' + str(tblCalc.scorePrint(genMotif)) + '.'
    #Onto next iteration.
    term += 1      
    
#Close data file.      
f.close()