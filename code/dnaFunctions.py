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


#==============================================================================================#
def allPossibleMotifs(length):

	motifs = map(''.join, itertools.product('ACGT', repeat=length))

	return motifs
