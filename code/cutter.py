import csv
import dnaFunctions
import sqlite3
import os
import ConfigParser

def main():

	parser = ConfigParser.SafeConfigParser()
	parser.read('config.ini')

	rootDirectory = parser.get('directories', 'root_dir')
	dataDirectory = os.path.join(rootDirectory, parser.get('directories', 'data_subdir'))
	bugDirectory = os.path.join(dataDirectory, parser.get('directories', 'bug_subdir'))
	databaseFileName = os.path.join(dataDirectory, parser.get('files', 'database_file_name'))

	dbConnection = sqlite3.connect(databaseFileName)
	cursor = dbConnection.cursor()

	sqlString = "SELECT name, type, subfolder, file_prefix FROM organisms WHERE type <> 'Human_Microbiom'"

	print '"organism_name","sequence_name","starting_location","ending_location","strand","type","sequence"'

	for row in cursor.execute(sqlString):

		name       = row[0].encode('ascii', 'ignore')
		category   = row[1].encode('ascii', 'ignore')
		subfolder  = row[2].encode('ascii', 'ignore')
		filePrefix = row[3].encode('ascii', 'ignore')

		fnaFile = filePrefix + ".fna.oneline"
		rntFile = filePrefix + ".rnt"
		pttFile = filePrefix + ".ptt"

		genomeFileName = os.path.join(bugDirectory, category, subfolder, fnaFile)
		rntFileName    = os.path.join(bugDirectory, category, subfolder, rntFile)
		pttFileName    = os.path.join(bugDirectory, category, subfolder, pttFile)

		data = generateSequenceMetaData(genomeFileName, rntFileName, pttFileName)

		for value in data.values():

			print '"' + name + '","' + value['name'] + '",' + str(value['start']) + ',' + str(value['end']) + ',"' + value['strand'] + '","' + value['type'] + '","' + value['sequence'] + '"'

	dbConnection.close()

def generateSequenceMetaData(genomeFileName, rntFileName, pttFileName):

		sequences = {}

		genomeFile = open(genomeFileName, 'r')
		rntFile = open(rntFileName, 'r')
		pttFile = open(pttFileName, 'r')

		rntMap = {'name':'Product', 'strand':'Strand', 'location':'Location'}
		dnaMap = {'name':'Gene', 'strand':'Strand', 'location':'Location'}

		readMetaDataFromFile(genomeFile, rntFile, 2, rntMap, 'RNA', sequences)
		readMetaDataFromFile(genomeFile, pttFile, 2, dnaMap, 'Gene', sequences)

		pttFile.close()
		rntFile.close()

		previousSequence = {'start': 0, 'end': 0, 'type': 'junk'}

		for startingLocation in sorted(sequences.iterkeys()):

			sequence = sequences[startingLocation]

			start  = int(previousSequence['end']) + 1
			end    = int(sequence['start']) - 1
			strand = '+'
			type   = 'DNA'
			name   = ''
			
			if start < end:

				sequenceStr = readFromFile(genomeFile, start, end, strand)
				dnaList  = {'start': start, 'end': end, 'strand': strand, 'type': type, 'name': name, 'sequence':sequenceStr}
				sequences[int(start)] = dnaList
				
			previousSequence = sequence

		return sequences


def readMetaDataFromFile(genomeFile, f, skipRows, columnMap, type, outputDictionary):

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
			name = ''

		sequence = readFromFile(genomeFile, begin, end, strand)

		outputDictionary[begin] = {'start':begin, 'end':end, 'strand':strand, 'type':type, 'name':name, 'sequence':sequence}


def readFromFile(readFile, start, end, strand):

	readLength = end - start + 1

	readFile.seek(start)
	sequence = readFile.read(readLength)

	if strand == '-':
		sequence = dnaFunctions.reverseComplement(sequence)

	return sequence


if __name__ == '__main__':
	main()