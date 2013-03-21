import dnaFunctions as DNA
import random, re

class Sequence:


	#==============================================================================================#
	def __init__(self, sequence, start=None, end=None, strand=None, type=None, name=None):

		self.sequence = sequence.strip()
		self.start    = start
		self.end      = end
		self.strand   = strand
		self.type     = type
		self.name     = name


	#==============================================================================================#
	def __str__(self):

		return self.sequence


	#==============================================================================================#
	def __len__(self):

		return len(self.sequence)


	#==============================================================================================#
	def __getitem__(self, key):

		section = 'X' + self.sequence
		return Sequence(section[key])


	#==============================================================================================#
	def reverseComplement(self):

		rComplement = DNA.reverseComplement(self.sequence)
		rComplement = Sequence(rComplement)

		return rComplement


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
				motifCounts[motif] = float(motifCounts[motif]) / len(self.sequence)

		return motifCounts


	#==============================================================================================#
	def getInvertedRepeats(self, min, max, mismatches):

		irs = {}

		motifs = self.motifs(min, max)

		# loop over each motif in the list of motifs
		for motif in motifs:

			# make a set to store the potential inverted repeats
			# using a set to avoid adding the same element multiple times
			potentialInvertedRepeats = set()

			# take the reverse complement of the motif
			reverseComplementMismatches = DNA.reverseComplementMismatches(motif, mismatches)

			# for each reverse complement with one or more mismatched base pairs
			for reverseComplement in reverseComplementMismatches:

				# if the reverse complement is also in the dictionary, then we know we have a possible IR
				if motifs.has_key(reverseComplement):

					# get the location(s) where this motif occurs
					motifLocations = motifs[motif]

					# get the location(s) where the reverse complement occurs
					reverseComplementLocations = motifs[reverseComplement]

					# loop over each starting position in the list of location(s) where the motif occurs
					for motifLocation in motifLocations:

						# loop over each starting position in the list of location(s) where the reverse complement occurs
						for reverseComplementLocation in reverseComplementLocations:

							# make sure there are no letters shared between the motifs
							if abs(reverseComplementLocation - motifLocation) >= len(motif):

								# construct a unique identifier for each motif and its reverse complement along with the starting locations
								# this is to avoid counting the same thing multiple times
								if motifLocation < reverseComplementLocation:
									identifier = motif + '_' + reverseComplement + '_' + str(motifLocation) + '_' + str(reverseComplementLocation)
								elif reverseComplementLocation < motifLocation:
									identifier = reverseComplement + '_' + motif + '_' + str(reverseComplementLocation) + '_' + str(motifLocation)
								else:
									identifier = None

								# add the identifier to the set of potential inverted repeats
								# the data structure is a set so if the identifer is already in the set, it will not be added again
								if identifier is not None:
									potentialInvertedRepeats.add(identifier)

			# get the length of the motif
			key = len(motif)

			# if we've already begun counting IRs of length 2 * key, then union the exisiting set with the generated set of potential inverted repeats
			if irs.has_key(key):
				irs[key] = irs[key].union(potentialInvertedRepeats)
			# otherwise, initialize an empty set
			else:
				irs[key] = set()

		return irs


	def getIRDistributions(self, min, max, mismatches, relativePercentages=True):

		irIdentifiers = self.getInvertedRepeats(min, max, mismatches)
		irDistributions = {}

		for key in irIdentifiers.keys():

			irSet = irIdentifiers[key]

			irDistribution = {}

			for invertedRepeat in irSet:

				regexString = '([ACGT]{' + str(key) + '})_[ACGT]{' + str(key) + '}_\d+_\d+'
				match = re.match(regexString, invertedRepeat)

				distributionKey = match.groups()[0]

				if irDistribution.has_key(distributionKey):
					irDistribution[distributionKey] = irDistribution[distributionKey] + 1
				else:
					irDistribution[distributionKey] = 1

			if relativePercentages == True:

				for irKey in irDistribution.keys():
					irDistribution[irKey] = float(irDistribution[irKey]) / len(self.sequence)

			irDistributions[key] = irDistribution

		return irDistributions


	def countInvertedRepeats(self, min, max, mismatches, relativePercentages=True):

		irIdentifiers = self.getInvertedRepeats(min, max, mismatches)

		countsOfIRs = {}

		for key in irIdentifiers.keys():

			countsOfIRs[key] = len(irIdentifiers[key])

			if relativePercentages == True:
				countsOfIRs[key] = float(countsOfIRs[key]) / len(self.sequence)

		for i in range(min, max+1):
			if not countsOfIRs.has_key(i):
				countsOfIRs[i] = 0

		return countsOfIRs


	#==============================================================================================#
	def shuffle(self):

		listOfChars = list(self.sequence)
		random.seed()

		for i in range(len(listOfChars) - 1, 0, -1):

			j = random.randint(0, i)

			temp = listOfChars[j]
			listOfChars[j] = listOfChars[i]
			listOfChars[i] = temp

		shuffledSequence = ''.join(listOfChars)

		return Sequence(shuffledSequence)
