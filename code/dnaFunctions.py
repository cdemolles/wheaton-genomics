import itertools

#==============================================================================================#
def complementBase(base):
	
		if base == 'A':
			return 'T'
		elif base == 'T':
			return 'A'
		elif base == 'G':
			return 'C'
		elif base == 'C':
			return 'G'
		else:
			return ''


#==============================================================================================#
def reverseComplement(motif):

	# initalize the reverse complement string
	reverseComplement = ''

	# for each base pair in the motif
	for base in motif:

		# find the complement of the base
		complement = complementBase(base)

		# add it to the beginning of the reverseComplement
		reverseComplement = complement + reverseComplement

	return reverseComplement


def reverseComplementMismatches(motif, numberOfMismatches):

	rComplement = reverseComplement(motif)

	if numberOfMismatches == 1:

		mismatches = set()

		for i in range(len(rComplement)):

			for basePair in allPossibleMotifs(1):

				reverseComplementList = list(rComplement)
				reverseComplementList[i] = basePair

				mismatch = ''.join(reverseComplementList)
				mismatches.add(mismatch)

		mismatches = list(mismatches)

		return mismatches

	else:

		return [rComplement]


#==============================================================================================#
def allPossibleMotifs(length):

	motifs = map(''.join, itertools.product('ACGT', repeat=length))

	return motifs
