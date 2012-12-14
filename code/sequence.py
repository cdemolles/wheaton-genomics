import dnaFunctions as DNA
import random

class Sequence:


	#==============================================================================================#
	def __init__(self, sequence, start=None, end=None, strand=None):

		self.sequence = sequence.strip()
		self.start    = start
		self.end      = end
		self.strand   = strand


	#==============================================================================================#
	def __str__(self):

		return self.sequence


	#==============================================================================================#
	def motifs(self, min, max):

		# tack on a dummy character to the sequence so the actual first character in the sequence
		# starts at index 1
		section = 'X' + self.sequence

		# start at the first character in the sequence
		start  = 1

		# initialize a blank dictionary to store all motifs of length min up to length max
		motifs = {}
	
		while start < len(section):
		
			strandLength = len(section)

			for i in range(min, max + 1): 

				if start <= strandLength - i:

					motif = section[start:start + i]

					if motifs.has_key(motif):
						motifs[motif].append(start)
					else:
						motifs[motif] = [start]
		
			start = start + 1
	
		return motifs


	#==============================================================================================#
	def countMotifs(self, min, max, relativePercentages=True):

		motifCounts = {}

		motifs = self.motifs(min, max)

		for motif in motifs.keys():

			motifCounts[motif] = len(motifs[motif])

			if relativePercentages == True:
				motifCounts[motif] = float(motifCounts[motif]) / (len(self.sequence) - len(motif) + 1)

		return motifCounts


	#==============================================================================================#
	def countInvertedRepeats(self, min, max, relativePercentages=True):

		countsOfIRs = {}

		motifs = self.motifs(min, max)

		for motif in motifs:

			count = 0

			# take the reverse complement of the motif
			reverseComplement = DNA.reverseComplement(motif)

			# if the reverse complement is also in the dictionary, then we know we have a possible IR
			if motifs.has_key(reverseComplement):

				motifLocations = motifs[motif]
				reverseComplementLocations = motifs[reverseComplement]

				for motifLocation in motifLocations:

					for reverseComplementLocation in reverseComplementLocations:

						if motifLocation < reverseComplementLocation:

							count = count + 1

			# get the length of the motif
			key = len(motif)

			if countsOfIRs.has_key(key):
				countsOfIRs[key] = countsOfIRs[key] + count
			else:
				countsOfIRs[key] = count

		for key in countsOfIRs.keys():

			countsOfIRs[key] = countsOfIRs[key] / 2

			if relativePercentages == True:
				countsOfIRs[key] = float(countsOfIRs[key]) / (len(self.sequence) - key + 1)

		for i in range(min, max+1):
			if not countsOfIRs.has_key(i):
				countsOfIRs[i] = 0

		return countsOfIRs


	def shuffle(self):

		listOfChars = list(self.sequence)
		random.shuffle(listOfChars)

		shuffledSequence = ''.join(listOfChars)

		return Sequence(shuffledSequence)
