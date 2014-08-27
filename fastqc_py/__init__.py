import sys
import numpy as np
import matplotlib.pyplot as plt

class fqFile:
	def __init__(self, inFile):
		"""
		Loads a fastq file for analysis
		
		Parameters:
		  inFile - File to be read. Only accepts fastq and fq files.
		"""
		extension = inFile.split('.')[-1]
		formats = ('fq','fastq')
		if extension in formats:
			self.inFile = inFile
		else:
			print "Bad file format: "+extension
			print "Requires either an fq or fastq file"
	def readLength(self, plot=False, printOut=True):
		"""
		Analyzes the fastq file for the sequence length

		Parameters:
		  plot     - produces a histogram plot of read lengths
		             [True, False] Default:False
		  printOut - prints statistics
		             [True, False] Default:True
		"""
		if not self.inFile:
			print "Needs an input file"
			return
		lens = []
		count = 0
		for seq,qual in fileGen(self.inFile):
			lens.append(len(seq))
			count += 1
		maxLen = max(lens)
		minLen = min(lens)
		if printOut:
			print "Finished reading %i reads" % (count)
			print "Min read length: %i" % (minLen)
			print "Max read length: %i" % (maxLen)
		if plot:
			plt.figure()
			plt.hist(lens)
			plt.title("Histogram of Read Lengths")
			plt.ylabel("Counts")
			plt.xlabel("Read Length")
			plt.show()
		self.maxLen = maxLen
		self.numReads = count

	def plotQual(self, printOut=True):
		"""
		View a boxplot of the qualities by base.

		Parameters:
		  printOut - print quality format
		             [True, False] Default:True
		"""
		if not 'self.maxLen' in locals():
			self.readLength(printOut=False)
		quals = initMatrix(self.maxLen)
		for seq,qual in fileGen(self.inFile):
			tmp = map(ord,qual)
			for i in xrange(len(tmp)):
				quals[i].append(tmp[i])
		qualRange = calcQualRange(quals)
		self.qualRange = qualRange
		print "Quality format: +%d" % (qualRange[0])
		plt.figure(figsize=(18,3))
		# switch to axes.bxp
		plt.boxplot(quals, sym='')
		plt.plot(range(1,len(quals)+1),map(np.mean, quals))
		plt.title("%s Quality Plot" % (self.inFile.split('/')[-1]))
		plt.ylabel("Quality Score")
		plt.ylim(qualRange)
		plt.tick_params(axis='x', which='both', labelbottom='off')
		plt.tight_layout()
		plt.show()

	def plotBaseBias(self):
		"""
		Plot the base bias by position.
		"""
		if not 'self.maxLen' in locals():
			self.readLength(printOut=False)
		bases = np.zeros((self.maxLen, 4))
		lookup = {'A':0,'G':1,'C':2,'T':3}
		nt = ('A','G','C','T')
		count = 0
		import time
		start = time.clock()
		for seq,qual in fileGen(self.inFile):
			for i in xrange(4):
				bases[seq==nt[i],i] += 1
			count += 1
			if not count % 100000: print "Finished %d of %d reads" % (count, self.numReads)
		print time.clock()-start
		sums = np.sum(bases,axis=1)
		plt.figure(figsize=(16,3))
		plt.plot(bases/np.matrix(sums).T)
		plt.legend(['A','G','C','T'],loc=5,bbox_to_anchor=(1.075,0.5))
		plt.title("%s Base Bias" % (self.inFile.split('/')[-1]))
		plt.ylabel("% of Bases")
		plt.subplots_adjust(left=0.035,right=0.93)
		plt.show()

def initMatrix(size):
	A = []
	for i in xrange(size):
		A.append([])
	return A

def fileGen(inFile):
	IF = open(inFile,'r')
	name1 = IF.readline()
	while name1:
		seq = np.core.defchararray.asarray(IF.readline().rstrip('\n'), itemsize=1)
		name1 = IF.readline()
		qual = IF.readline().rstrip('\n')
		name1 = IF.readline()
		yield((seq,qual))
	IF.close()

def calcQualRange(quals):
	minQual = min(map(min, quals))
	maxQual = max(map(max, quals))
	if maxQual > 74:
		return (64, 104)
	else:
		return (33, 74)
