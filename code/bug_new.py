import os
import csv
import glob
import dnaFunctions
from sequence import *

class Bug:


	#==============================================================================================#
	def __init__(self, name):

		rootDirectory = '/home/chris/Documents/microbes/HMP_new/'

		self.name = name
		self.rna  = None
		self.dna  = None
		self.sequenceMetaData = None

		self.genomeFileName   = glob.glob(rootDirectory + name + '*/*_complete_genome.fna.oneline')[0]
		self.pttFileName      = glob.glob(rootDirectory + name + '*/*_complete_genome.ptt')[0]
		self.rntFileName      = glob.glob(rootDirectory + name + '*/*_complete_genome.rnt')[0]

		f = open(self.genomeFileName)
		genomeStr = f.readline()
		genomeStr = genomeStr.strip()
		genomeStr = genomeStr[1:]
		self.genome = Sequence(genomeStr)


	def __getitem__(self, key):

		sequenceSlice = self.genome[key]
		return sequenceSlice


	#==============================================================================================#
	# Returns a dictionary of sequences
	#==============================================================================================#
	def getRna(self):

		if self.rna is None:

			if self.sequenceMetaData is None:
				self.generateSequenceMetaData()

			self.rna = self.cacheSequences('RNA')

		return self.rna


	#==============================================================================================#
	def getDna(self):

		if self.dna is None:
			
			if self.sequenceMetaData is None:
				self.generateSequenceMetaData()

			self.dna = self.cacheSequences('DNA')

		return self.dna


	#==============================================================================================#
	def shuffle(self):

		shuffledBug = Bug(self.name)
		shuffledBug.genome = shuffledBug.genome.shuffle()

		return shuffledBug


	#==============================================================================================#
	def readFromFile(self, readFile, start, end, strand):

		readLength = end - start + 1

		readFile.seek(start)
		sequence = readFile.read(readLength)

		if strand == '-':
			sequence = dnaFunctions.reverseComplement(sequence)

		return sequence


	#==============================================================================================#
	def cacheSequences(self, type):

		listOfSequences = {}

		#genomeFile = open(self.genomeFileName, 'r')

		for value in self.sequenceMetaData.values():

			if value['type'] == type:

				sequence = self.genome[value['start']:value['end']+1]

				if value['strand'] == '-':
					sequence = sequence.reverseComplement()

				sequence.start = value['start']
				sequence.end   = value['end']
				sequence.strand = value['strand']

				if type == 'DNA':

					if not listOfSequences.has_key(value['start']):
						listOfSequences[value['start']] = sequence

				else:

					if not listOfSequences.has_key(value['name']):
						listOfSequences[value['name']] = sequence

		#genomeFile.close()

		return listOfSequences


	#==============================================================================================#
	def generateSequenceMetaData(self):

		sequences = {}

		rntFile = open(self.rntFileName, 'r')
		pttFile = open(self.pttFileName, 'r')

		rntMap = {'name':'Product', 'strand':'Strand', 'location':'Location'}
		dnaMap = {'name':'Gene', 'strand':'Strand', 'location':'Location'}

		self.readMetaDataFromFile(rntFile, 2, rntMap, 'RNA', sequences)
		self.readMetaDataFromFile(pttFile, 2, dnaMap, 'Gene', sequences)

		pttFile.close()
		rntFile.close()

		previousSequence = {'start': 0, 'end': 0, 'type': 'junk'}

		for startingLocation in sorted(sequences.iterkeys()):

			sequence = sequences[startingLocation]

			start  = int(previousSequence['end']) + 1
			end    = int(sequence['start']) - 1
			strand = sequence['strand']
			type   = 'DNA'
			name   = 'unknown'
			
			if start < end:

				dnaList  = {'start': start, 'end': end, 'strand': strand, 'type': type, 'name': name}

				sequences[int(start)] = dnaList
				previousSequence = sequence

		self.sequenceMetaData = sequences


	#==============================================================================================#
	def readMetaDataFromFile(self, f, skipRows, columnMap, type, outputDictionary):

		f.seek(0)

		for i in range(skipRows):
			f.readline()

		for line in csv.DictReader(f, delimiter='\t'):

			name   = line[columnMap['name']]
			strand = line[columnMap['strand']]
			begin  = int(line[columnMap['location']].split('..')[0])
			end    = int(line[columnMap['location']].split('..')[1])

			name = name.split(' ')[0]

			if name == '-':
				name = 'unknown'

			outputDictionary[begin] = {'start':begin, 'end':end, 'strand':strand, 'type':type, 'name':name}

