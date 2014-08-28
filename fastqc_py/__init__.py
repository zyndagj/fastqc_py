import sys
import numpy as np
import matplotlib.pyplot as plt
import time
from multiprocessing import Process, Pipe, cpu_count

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
		count = 0
		for seq,qual in fileGen(self.inFile):
			tmp = map(ord,qual)
			for i in xrange(len(tmp)):
				quals[i].append(tmp[i])
			count += 1
		print "Read %i reads" % (count)
		outMatrix = np.zeros((self.maxLen, 105), dtype=np.uint32)
		for i in xrange(len(quals)):
			for qual in quals[i]:
				outMatrix[i,qual] += 1
		return outMatrix
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

	def plotBXP(self, printOut=True):
		"""
		View a boxplot of the qualities by base.

		Parameters:
		  printOut - print quality format
		             [True, False] Default:True
		"""
		if not 'self.maxLen' in locals():
			self.readLength(printOut=False)
		qualMatrix = np.zeros((self.maxLen, 105), dtype=np.uint32)
		tmpMatrix = np.zeros((self.maxLen, 105), dtype=np.uint32)
		for seq,qual in fileGen(self.inFile):
			tmp = map(ord,qual)
			qualMatrix[(range(len(tmp)),tmp)] += 1
		return qualMatrix
		self.qualRange = calcQualRangeBXP(qualMatrix)
		print "Quality format: +%d" % (self.qualRange[0])
		sys.exit()
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

	def plotBaseBias(self, nCores=1, verbose=False):
		"""
		Plots the base bias, or per base sequence content.

		Parameters
		=============
		nCores	INT	Number of cpu cores to use. -1 to use all cores. (Default: 1)
		verbose	BOOL	Print runtime (Default: False)
		"""
		if not 'self.maxLen' in locals():
			self.readLength(printOut=False)
		bases = np.zeros((self.maxLen, 4))
		if nCores == 1:
			nt = ('A','G','C','T')
			count = 0
			start = time.clock()
			for seq,qual in fileGen(self.inFile):
				npSeq = np.core.defchararray.asarray(seq, itemsize=1)
				for i in xrange(4):
					bases[npSeq==nt[i],i] += 1
				count += 1
				if not count % 100000: print "Finished %d of %d reads" % (count, self.numReads)
			if verbose:
				print "CPU time: %.3f seconds" % (time.clock()-start)
		else:
			maxCores = cpu_count()
			if nCores > maxCores:
				print "Got %i cores, but system only has %i cores. Using %i cores." % (nCores, maxCores, maxCores)
				nCores = maxCores
			if nCores == -1:
				nCores = maxCores
			if verbose:
				print "Using %i cores" % (nCores)
			p = [] # process array
			pConns = [] # parent connection array
			for i in xrange(nCores):
				pConn, cConn = Pipe() #returns (parent connection, child connection)
				pConns.append(pConn)
				# initialize processes
				p.append(Process(target=bbWorker, args=(self.inFile, self.maxLen, i, nCores, cConn)))
			wallStart = time.time()
			for i in xrange(nCores): # start processes
				p[i].start()
			cpuTotal = 0
			for i in xrange(nCores):
				tmpBases, cpuTime = pConns[i].recv() # get results from processes
				bases += tmpBases
				cpuTotal += cpuTime
			wallTime = time.time()-wallStart
			for i in xrange(nCores):
				p[i].join()
			if verbose:
				print "CPU time: %.3f seconds" % (cpuTotal)
				print "Walltime: %.3f seconds" % (wallTime)
		sums = np.sum(bases,axis=1)
		plt.figure(figsize=(12,4))
		plt.plot(bases/np.matrix(sums).T)
		plt.legend(['A','G','C','T'],loc=5,bbox_to_anchor=(1.1,0.5))
		plt.title("%s Base Bias" % (self.inFile.split('/')[-1]))
		plt.ylabel("% of Bases")
		plt.subplots_adjust(left=0.05,right=0.91)
		plt.show()

def bbWorker(inFile, maxLen, wid, procs, cconn):
	"""
	base bias worker called by plotBaseBias for parallel computation
	"""
	nt = ('A','G','C','T')
	tmpBases = np.zeros((maxLen,4))
	count = 0
	myCount = 0
	cpuStart = time.clock()
	for seq, qual in fileGen(inFile):
		if count % procs == wid:
			npSeq = np.core.defchararray.asarray(seq, itemsize=1)
			for i in xrange(4):
				tmpBases[npSeq==nt[i],i] += 1
			myCount += 1
		count += 1
	cpuTime = time.clock()-cpuStart
	cconn.send((tmpBases,cpuTime))
	cconn.close()

def initMatrix(size):
	A = []
	for i in xrange(size):
		A.append([])
	return A

def fileGen(inFile):
	IF = open(inFile,'r')
	name1 = IF.readline()
	while name1:
		seq = IF.readline().rstrip('\n')
		name1 = IF.readline()
		qual = IF.readline().rstrip('\n')
		name1 = IF.readline()
		yield((seq,qual))
	IF.close()

def calcQualRangeBXP(qualMatrix):
	base, qual = np.where(qualMatrix > 0)
	minQual = min(qual)
	maxQual = max(qual)
	if maxQual > 74:
		return (64, 104)
	else:
		return (33, 74)

def calcQualRange(quals):
	minQual = min(map(min, quals))
	maxQual = max(map(max, quals))
	if maxQual > 74:
		return (64, 104)
	else:
		return (33, 74)
