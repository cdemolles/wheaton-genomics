from bug import *
from sequence import *
import glob
import os
import re
import dnaFunctions

def main():

	svmOutputRNA = open('../data/SVM_RNA_mismatch.tsv', 'w')
	svmOutputDNA = open('../data/SVM_DNA_mismatch.tsv', 'w')
	
	allPossibleFourMers = dnaFunctions.allPossibleMotifs(4)

	outputHeading(svmOutputRNA, allPossibleFourMers)
	outputHeading(svmOutputDNA, allPossibleFourMers)

	for filePath in glob.glob('/home/chris/Documents/microbes/HMP_new/*'):

		 bugName = os.path.split(filePath)[1]
		 bugName = re.split('_uid', bugName)[0]

		 print bugName

		 if glob.glob('/home/chris/Documents/microbes/HMP_new/' + bugName + '*/*_complete_genome*') != []:

			 bug = Bug(bugName)

			 # gets a dictionary of RNA sequences (with the key as the name of the RNA and the value as a sequence object)
			 bugRNA = bug.getRna()

			 # gets a dictionary of DNA sequences (with the key as the starting location of the DNA and the value as a sequence object)
			 bugDNA = bug.getDna()

			 outputData(svmOutputRNA, bugRNA, bugName, 'RNA', allPossibleFourMers)
			 outputData(svmOutputDNA, bugDNA, bugName, 'DNA', allPossibleFourMers)

	svmOutputRNA.close()
	svmOutputDNA.close()


def outputHeading(outputFile, allPossibleFourMers):

	outputFile.write('Name')
	outputFile.write('\t')

	for motif in allPossibleFourMers:
		outputFile.write(motif)
		outputFile.write('\t')

	outputFile.write('IRs of Length 6')
	outputFile.write('\t')

	outputFile.write('IRs of Length 8')
	outputFile.write('\t')

	outputFile.write('IRs of Length 10')
	outputFile.write('\t')

	outputFile.write('IRs 1 Mismatch Length 6')
	outputFile.write('\t')

	outputFile.write('IRs 1 Mismatch Length 8')
	outputFile.write('\t')

	outputFile.write('Type')
	outputFile.write('\n')

def outputData(outputFile, data, bugName, type, allPossibleFourMers):

	for sequenceName, sequence in data.iteritems():

		motifs          = sequence.countMotifs(4, 4)
		invertedRepeats = sequence.countInvertedRepeats(3, 5, 0)
		invertedRepeatsMismatch = sequence.countInvertedRepeats(3, 4, 1)

		strandName = ''

		if sequence.strand == '+':
			strandName = 'direct'
		elif sequence.strand == '-':
			strandName = 'indirect'

		outputName = ''

		if type == 'DNA':
			outputName = 'XXX'
		else:
			outputName = sequenceName

		outputName = bugName + '_' + type + '_' + outputName + '_' + str(sequence.start) + '_' + str(sequence.end) + '_' + strandName

		outputFile.write(outputName)
		outputFile.write('\t')

		for motif in allPossibleFourMers:

			if motifs.has_key(motif):
				strCount = str(motifs[motif])
				outputFile.write(strCount)
			else:
				outputFile.write('0')

			outputFile.write('\t')

		outputFile.write(str(invertedRepeats[3]))
		outputFile.write('\t')
		outputFile.write(str(invertedRepeats[4]))
		outputFile.write('\t')
		outputFile.write(str(invertedRepeats[5]))
		outputFile.write('\t')
		outputFile.write(str(invertedRepeatsMismatch[3]))
		outputFile.write('\t')
		outputFile.write(str(invertedRepeatsMismatch[4]))
		outputFile.write('\t')

		if type == 'DNA':
			outputFile.write('notRNA')
		else:
			outputFile.write(type)

		outputFile.write('\n')

if __name__ == '__main__':
	main()
