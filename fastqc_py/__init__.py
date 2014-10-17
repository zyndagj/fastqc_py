import sys
import os.path
import numpy as np
import matplotlib.pyplot as plt
import time
from multiprocessing import Process, Pipe, cpu_count
from array import array
from collections import Counter
import matplotlib.cm as cm
from Bio import pairwise2, SeqIO, Seq
from cIO import cFileGen

baseDict = {'A':'0','G':'1','C':'2','T':'3'}
bases = ('A','G','C','T','N')

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
		#for seq,qual in fileGen(self.inFile):
		for seq,qual in cFileGen(self.inFile):
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

	def plotQual(self, printOut=True, nCores=1, verbose=False, plot=True):
		"""
		View a boxplot of the qualities by base.

		Parameters
		=======================
		printOut	BOOL	print quality format (Default: True)
		nCores		INT	Number of cpu cores to use. -1 to use all cores. (Default: 1)
		verbose		BOOL	Print the runtime (Default: False)
		"""
		if not 'self.maxLen' in locals():
			self.readLength(printOut=False)
		quals = initMatrix(self.maxLen)
		if nCores == 1:
			cpuStart = time.clock()
			wallStart = time.time()
			for seq,qual in fileGen(self.inFile):
				tmp = map(ord,qual)
				for i in xrange(len(tmp)):
					quals[i].append(tmp[i])
			cpuTotal = time.clock() - cpuStart
			wallTime = time.time() - wallStart
		else:
			nCores = setCores(nCores, verbose)
			p = [] # process array
			pConns = [] # parent connection array
			for i in xrange(nCores):
				pConn, cConn = Pipe() #returns (parent connection, child connection)
				pConns.append(pConn)
				p.append(Process(target=qualWorker, args=(self.inFile, self.maxLen, i, nCores, cConn))) # initialize processes
			wallStart = time.time()
			for i in xrange(nCores): # start processes
				p[i].start()
			cpuTotal = 0
			for i in xrange(nCores):
				tmpQuals, cpuTime = pConns[i].recv() # get results from processes
				for i in xrange(len(tmpQuals)):
					quals[i].extend(tmpQuals[i])
				cpuTotal += cpuTime
			wallTime = time.time()-wallStart
			for i in xrange(nCores):
				p[i].join()
			print sum(map(sum, quals))
		if verbose:
			print "CPU time: %.3f seconds" % (cpuTotal)
			print "Walltime: %.3f seconds" % (wallTime)
		qualRange = calcQualRange(quals)
		self.qualRange = qualRange
		if printOut:
			print "Quality format: +%d" % (qualRange[0])
		if plot:
			plt.figure(figsize=(18,3))
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
		bArray = np.zeros((self.maxLen,5), dtype=np.uint32)
		if nCores == 1:
			count = 0
			start = time.clock()
			#for seq,qual in fileGen(self.inFile):
			for seq,qual in cFileGen(self.inFile):
				npSeq = np.core.defchararray.asarray(seq, itemsize=1)
				for i in xrange(5):
					bArray[npSeq==bases[i],i] += 1
				count += 1
				if not count % 100000: print "Finished %d of %d reads" % (count, self.numReads)
			if verbose:
				print "CPU time: %.3f seconds" % (time.clock()-start)
		else:
			nCores = setCores(nCores, verbose)
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
				bArray += tmpBases
				cpuTotal += cpuTime
			wallTime = time.time()-wallStart
			for i in xrange(nCores):
				p[i].join()
			if verbose:
				print "CPU time: %.3f seconds" % (cpuTotal)
				print "Walltime: %.3f seconds" % (wallTime)
		sums = np.sum(bArray,axis=1)
		plt.figure(figsize=(12,4))
		for i in range(5):
			plt.plot(bArray[:,i]/sums)
		plt.legend(bases,loc=5,bbox_to_anchor=(1.1,0.5))
		plt.title("%s Base Bias" % (self.inFile.split('/')[-1]))
		plt.ylabel("% of Bases")
		plt.subplots_adjust(left=0.05,right=0.91)
		plt.show()

	def calcKmers(self, k=5, nCores=1, plot=True, adapt=True, verbose=False):
		"""
		Calculate and make a scatter plot of k-mers in reads.

		Parameters
		======================
		k	INT	k-mer size (Default: 5)
		nCores	INT	cores to use (Default: 1)
		plot	BOOL	plot the output (Default: True)
		adapt	BOOL	align kmers against illumina adapters (Default:True)
		verbose	BOOL	print the runtime (Default: False)
		"""
		kmerDict = Counter()
		if nCores == 1:
			wallStart = time.time()
			cpuStart = time.clock()
			for seq, qual in fileGen(self.inFile):
				tmp = [seq[i:i+k] for i in xrange(0,len(seq)-(k-1),3)]
				for i in tmp:
					if not 'N' in i:
						kmerDict[i]+= 1
			wallTime = time.time()-wallStart
			cpuTotal = time.clock()
		else:
			nCores = setCores(nCores, verbose)
			p = [] # process array
			pConns = [] # parent connection array
			for i in xrange(nCores):
				pConn, cConn = Pipe() #returns (parent connection, child connection)
				pConns.append(pConn)
				# initialize processes
				p.append(Process(target=kmerWorker, args=(self.inFile, k, i, nCores, cConn)))
			wallStart = time.time()
			for i in xrange(nCores): # start processes
				p[i].start()
			cpuTotal = 0
			for i in xrange(nCores):
				kmerTop100, cpuTime = pConns[i].recv() # get results from processes
				cpuTotal += cpuTime
				for k,v in kmerTop100:
					kmerDict[k]+=v
			wallTime = time.time()-wallStart
			for i in xrange(nCores):
				p[i].join()
		if verbose:
			print "CPU time: %.3f seconds" % (cpuTotal)
			print "Walltime: %.3f seconds" % (wallTime)
		top20 = kmerDict.most_common(20)
		vals = map(lambda y: y[1], top20)
		bottoms = np.cumsum([0]+vals[:-1])
		if plot:
				plt.figure(figsize=(4,8))
				plt.axis('off')
				plt.bar(np.zeros(20), vals, width=np.ones(20), color=cm.Set1(np.linspace(0,1,20)), bottom=bottoms)
				plt.xlim((-0.05,1.6))
				for i in xrange(20):
					plt.text(1.1, bottoms[i]+vals[i]/2.0, top20[i][0], verticalalignment='center', family='monospace')
					plt.text(0.5, bottoms[i]+vals[i]/2.0, str(vals[i]), va='center',ha='center')
				fName = self.inFile.split('/')[-1]
				plt.title("Top 20 K-mers in "+fName)
				plt.tight_layout()
				plt.show()
		aFile = os.path.join(os.path.dirname(__file__),'adapter_sequences.fa')
		if adapt:
			aCounter = Counter()
			IF = open(aFile,'r')
			for record in SeqIO.parse(IF,'fasta'):
				for i in range(20):
					results = pairwise2.align.localms(top20[i][0],str(record.seq),2,-1,-2.0,-0.1)
					if results:
						if results[0][2] > 16:
							aCounter[record.name]+=1
			print("%-35s %-10s %s" % ("Adapter","Num Hits","Sequence"))
			rec_dict = SeqIO.index(aFile,'fasta')
			for k,v in aCounter.most_common(10):
				
				print("%-35s %-10d %s"%(k,v,str(rec_dict[k].seq)))

#	def plotBXP(self, printOut=True, verbose=False):
#		"""
#		View a boxplot of the qualities by base.
#
#		Parameters
#		=======================
#		printOut	BOOL	 print quality format (Default: True)
#		verbose		BOOL	Print the runtime (Default: False)
#		"""
#		if not 'self.maxLen' in locals():
#			self.readLength(printOut=False)
#		qualArray = np.zeros(self.maxLen*105, dtype=np.uint32)
#		start = time.clock()
#		for seq,qual in fileGen(self.inFile):
#			tmp = map(ord,qual)
#			qualArray[np.arange(len(tmp))*105+tmp] += 1
#		if verbose:
#			print "CPU time: %.3f seconds" % (time.clock()-start)
#		qualMatrix = qualArray.reshape((self.maxLen,105))
#		self.qualRange, boxList = processMatrixBXP(qualMatrix)
#		print "Quality format: +%d" % (self.qualRange[0])
#		plt.figure(figsize=(18,3))
#		fig, ax = plt.subplots(111)
#		ax[0].bxp(boxList)
#		# switch to axes.bxp
#		#plt.boxplot(quals, sym='')
#		#plt.plot(range(1,len(quals)+1),map(np.mean, quals))
#		plt.title("%s Quality Plot" % (self.inFile.split('/')[-1]))
#		plt.ylabel("Quality Score")
#		plt.ylim(qualRange)
#		plt.tick_params(axis='x', which='both', labelbottom='off')
#		plt.tight_layout()
#		plt.show()

def bbWorker(inFile, maxLen, wid, procs, cconn):
	"""
	base bias worker called by plotBaseBias for parallel computation
	"""
	tmpBases = np.zeros((maxLen,5), dtype=np.uint32)
	count = 0
	myCount = 0
	cpuStart = time.clock()
	#for seq, qual in fileGen(inFile):
	for seq, qual in cFileGen(inFile):
		if count % procs == wid:
			npSeq = np.core.defchararray.asarray(seq, itemsize=1)
			for i in xrange(5):
				tmpBases[npSeq==bases[i],i] += 1
			myCount += 1
		count += 1
	cpuTime = time.clock()-cpuStart
	cconn.send((tmpBases,cpuTime))
	cconn.close()

def qualWorker(inFile, maxLen, wid, procs, cConn):
	"""
	Quality worker called by plotQuality
	"""
	count = 0
	quals = initMatrix(maxLen)
	cpuStart = time.clock()
	for seq,qual in fileGen(inFile):
		if count % procs == wid:
			tmp = map(ord,qual)
			for i in xrange(len(tmp)):
				quals[i].append(tmp[i])
		count += 1
	cpuTime = time.clock()-cpuStart
	cConn.send((quals, cpuTime))
	cConn.close()

def kmerWorker(inFile,k,wid,procs,cConn):
	kmerDict = Counter()
	cpuStart = time.clock()
	count = 0
	for seq, qual in fileGen(inFile):
		if count % procs == wid:
			#tmp = [seq[i:i+k] for i in xrange(len(seq)-(k-1))]
			tmp = [seq[i:i+k] for i in xrange(0,len(seq)-(k-1),3)]
			for i in tmp:
				if not 'N' in i:
					kmerDict[i]+= 1
		count += 1
	cpuTime = time.clock()
	cConn.send((kmerDict.most_common(100),cpuTime))
	cConn.close()

def initMatrix(size):
	A = []
	for i in xrange(size):
		A.append(array('B',[]))
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

def processMatrixBXP(qualMatrix):
	base, qual = np.where(qualMatrix > 0)
	minQual = min(qual)
	maxQual = max(qual)
	if maxQual > 74:
		qualRange = (64, 104)
	else:
		qualRange = (33, 74)
	boxDictList = map(lambda y: calcBox(y,qualRange),qualMatrix)
	return qualRange, boxDictList

def calcBox(array, qualRange):
	numVals = sum(array)
	q1 = reducedPercentile(array, 25)
	q2 = reducedPercentile(array, 50)
	q3 = reducedPercentile(array, 75)
	IQR = q3-q1
	wLow = q1-1.5*IQR
	if wLow < qualRange[0]:
		wLow = qualRange[0]
	wHigh = q3+1.5*IQR
	if wHigh > qualRange[1]:
		wHigh = qualRange[1]
	return {'q1':q1, 'med':q2, 'q3':q3, 'whislo':wLow, 'whishi':wHigh}

def reducedPercentile(array, percent):
	location = sum(array)*float(percent)/100.0
	accum = 0
	if location != int(location):
		for i in xrange(len(array)):
			accum += array[i]
			if accum > location:
				return i
	else:
		for i in xrange(len(array)):
			for k in xrange(array[i]):
				accum += 1
				if accum == location:
					first = i
				elif accum == location+1:
					second = i
					return (first+second)/2.0

def calcQualRange(quals):
	minQual = min(map(min, quals))
	maxQual = max(map(max, quals))
	if maxQual > 74:
		return (64, 104)
	else:
		return (33, 74)

def setCores(nCores, verbose):
	tmpCores = nCores
	maxCores = cpu_count()
	if nCores > maxCores:
		print "Got %i cores, but system only has %i cores. Using %i cores." % (nCores, maxCores, maxCores)
		tmpCores = maxCores
	elif nCores == -1:
		tmpCores = maxCores
	if verbose:
		print "Using %i cores" % (tmpCores)
	return tmpCores
