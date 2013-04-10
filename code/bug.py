import os, glob, dnaFunctions, sqlite3
from sequence import Sequence
from ConfigParser import SafeConfigParser


class Bug:


	#-----------------------------------------------------------------------------------------------------------------#
	#==============================================================================================#
	# Constructor
	# SUMMARY: Initializes bug object
	# PRE: assigned(name)
	# POST: Initializes a new bug object with name and genome based on that name
	#==============================================================================================#
	def __init__(self, name):

		# use a parser to open the configuration file
		parser = SafeConfigParser()
		parser.read('config.ini')

		# get the data directory from the configuration file
		dataDirectory = os.path.join(parser.get('directories', 'root_dir'), parser.get('directories', 'data_subdir'))

		# get the bug directory from the configuration file
		bugDirectory  = os.path.join(dataDirectory, parser.get('directories', 'bug_subdir'))

		# get the database file name and the name of the bug
		self.databaseFileName = os.path.join(dataDirectory, parser.get('files', 'database_file_name'))
		
		dbConnection = sqlite3.connect(self.databaseFileName)
		cursor = dbConnection.cursor()
		parameters = (name,)

		bugDataQuery = "SELECT name, type, subfolder, file_prefix, kingdom, cellular, category FROM organisms WHERE name = ?"

		cursor.execute(bugDataQuery, parameters)

		bugData = cursor.fetchone()

		if bugData is not None:

			genomeFileName = os.path.join(bugDirectory + '/' + bugData[1].encode('ascii', 'ignore') + '/' + bugData[2].encode('ascii', 'ignore') + '/' + bugData[3].encode('ascii', 'ignore') + '.fna.oneline')
			
			self.bugName  = bugData[0].encode('ascii', 'ignore')
			self.kingdom  = bugData[4]#.encode('ascii', 'ignore')
			self.cellular = bugData[5]#.encode('ascii', 'ignore')
			self.group    = bugData[6]#.encode('ascii', 'ignore')
			self.genome   = Sequence(self.loadGenome(genomeFileName))
			
		else:
			print 'Invalid bug name. Check to make sure the bug name is spelled correctly or its genome file exists in its directory.'
		
	#-----------------------------------------------------------------------------------------------------------------#


	#-----------------------------------------------------------------------------------------------------------------#
	#==============================================================================================#
	# __getitem__()
	# SUMMARY: Overloads indexing and slicing operator
	# PRE: assigned(self.genome), assigned(key)
	# POST: Returns self.genome[key]
	#==============================================================================================#
	def __getitem__(self, key):

		# indexes or slices the genome and returns it
		sequenceSlice = self.genome[key]
		return sequenceSlice
	#-----------------------------------------------------------------------------------------------------------------#


	#-----------------------------------------------------------------------------------------------------------------#
	#==============================================================================================#
	# __str__()
	# SUMMARY: Overloads printing function
	# PRE: assigned(self.genome)
	# POST: Returns the genome as a string
	#==============================================================================================#
	def __str__(self):

		return self.genome.sequence
	#-----------------------------------------------------------------------------------------------------------------#


	#-----------------------------------------------------------------------------------------------------------------#
	#==============================================================================================#
	# __len__()
	# SUMMARY: Returns the length of the bug's genome
	# PRE: assigned(self.genome)
	# POST: Returns the length of self.genome
	#==============================================================================================#
	def __len__(self):

		return len(self.genome)
	#-----------------------------------------------------------------------------------------------------------------#


	#-----------------------------------------------------------------------------------------------------------------#
	#==============================================================================================#
	# getSequences()
	# SUMMARY: Retrieves a dictionary of sequences of a user-defined type
	# PRE: assigned(self.bugName), assigned(self.databaseName), assigned(type)
	# POST: Returns a dictionary with the starting position of the sequence as the key and the sequence as the value
	#==============================================================================================#
	def getSequences(self, type):

		# initialize an empty dictionary to hold the sequences
		sequences = {}

		# connect to the database
		dbConnection = sqlite3.connect(self.databaseFileName)

		# make a database cursor
		cursor = dbConnection.cursor()

		# store the parameters for the query
		parameters = (self.bugName, type,)

		# extract the sequence data from the database given the name of the bug and the type of sequence
		sqlString = "SELECT sequence, starting_location, ending_location, strand, type, sequence_name FROM sequences WHERE organism_name = ? AND type = ?"

		for row in cursor.execute(sqlString, parameters):

			sequenceStr = row[0].encode('ascii', 'ignore')
			start       = row[1]
			end         = row[2]
			strand      = row[3].encode('ascii', 'ignore')
			type        = row[4].encode('ascii', 'ignore')
			name        = row[5]

			if name is None or len(name) == 0:
				name = 'unknown'
			else:
				name = name.encode('ascii', 'ignore')
				
			sequence = Sequence(sequenceStr, start, end, strand, type, name)
			sequences[start] = sequence

		dbConnection.close()

		return sequences
	#-----------------------------------------------------------------------------------------------------------------#


	#==============================================================================================#
	# Returns a dictionary of sequences
	#==============================================================================================#
	def getUniqueRNA(self):
		
		sequences = {}

		dbConnection = sqlite3.connect(self.databaseFileName)
		cursor = dbConnection.cursor()

		parameters = (self.bugName, 'RNA',)

		sqlString = "SELECT sequence, starting_location, ending_location, strand, type, sequence_name FROM sequences WHERE organism_name = ? AND type = ? GROUP BY sequence_name"

		for row in cursor.execute(sqlString, parameters):

			sequenceStr = row[0].encode('ascii', 'ignore')
			start       = row[1]
			end         = row[2]
			strand      = row[3].encode('ascii', 'ignore')
			type        = row[4].encode('ascii', 'ignore')
			name        = row[5]

			if name is None or name == '':
				name = 'unknown'
			else:
				name = name.encode('ascii', 'ignore')
				
			sequence = Sequence(sequenceStr, start, end, strand, type, name)
			sequences[name] = sequence

		dbConnection.close()

		return sequences


	#==============================================================================================#
	def getRNA(self):

		return self.getSequences('RNA')


	#==============================================================================================#
	def getDNA(self):

		return self.getSequences('DNA')


	#==============================================================================================#
	def getAnnotation(self, start, end):

		annotationList = []

		dbConnection = sqlite3.connect(self.databaseFileName)
		cursor = dbConnection.cursor()

		parameters = (self.bugName, start, end, start, end, start, end)

		for row in cursor.execute("SELECT sequence, starting_location, ending_location, strand, type, sequence_name FROM sequences WHERE organism_name = ? AND (((? BETWEEN starting_location AND ending_location) OR (? BETWEEN starting_location AND ending_location)) OR ((starting_location BETWEEN ? AND ?) OR (ending_location BETWEEN ? AND ?))) ORDER BY starting_location, ending_location", parameters):

			sequenceStr = row[0].encode('ascii', 'ignore')
			start       = row[1]
			end         = row[2]
			strand      = row[3].encode('ascii', 'ignore')
			type        = row[4].encode('ascii', 'ignore')
			name        = row[5]

			if name is None:
				name = 'unknown'
			else:
				name = name.encode('ascii', 'ignore')
				
			sequence = Sequence(sequenceStr, start, end, strand, type, name)

			annotationList.append(sequence)

		dbConnection.close()

		return annotationList


	#==============================================================================================#
	def shuffle(self):

		shuffledBug = Bug(self.bugName)
		shuffledBug.genome = shuffledBug.genome.shuffle()

		return shuffledBug


	#==============================================================================================#
	def printReport(self):

		rna = self.getRna()
		dna = self.getDna()
		
		print '=================================================================================='
		print 'The name of this bug is', self.name
		print 'The length of its genome is', len(self.genome)
		print ''

		for rnaName, sequence in rna.iteritems():

			perfectIRs = sequence.countInvertedRepeats(3, 4, 0)
			IRsOneMismatch = sequence.countInvertedRepeats(3, 4, 1)

			if sequence.strand == '+':
				strand = 'direct'
			else:
				strand = 'indirect'

			print 'There is a', rnaName, 'RNA on the', strand, 'strand that begins at', sequence.start, 'and ends at', sequence.end
			print 'The length of the sequence is', len(sequence)
			print 'The ratio of potential perfect IRs of length 6 to the length of the sequence is', perfectIRs[3]
			print 'The ratio of potential perfect IRs of length 8 to the length of the sequence is', perfectIRs[4]
			print 'The ratio of potential IRs of length 6 with 1 mismatch to the length of the sequence is', IRsOneMismatch[3]
			print 'The ratio of potential IRs of length 8 with 1 mismatch to the length of the sequence is', IRsOneMismatch[4]
			print ''

		for dnaStart, sequence in dna.iteritems():

			perfectIRs = sequence.countInvertedRepeats(3, 4, 0)
			IRsOneMismatch = sequence.countInvertedRepeats(3, 4, 1)

			if sequence.strand == '+':
				strand = 'direct'
			else:
				strand = 'indirect'

			print 'There is a DNA sequence on the', strand, 'strand that begins at', sequence.start, 'and ends at', sequence.end
			print 'The length of the sequence is', len(sequence)
			print 'The ratio of potential perfect IRs of length 6 to the length of the sequence is', perfectIRs[3]
			print 'The ratio of potential perfect IRs of length 8 to the length of the sequence is', perfectIRs[4]
			print 'The ratio of potential IRs of length 6 with 1 mismatch to the length of the sequence is', IRsOneMismatch[3]
			print 'The ratio of potential IRs of length 8 with 1 mismatch to the length of the sequence is', IRsOneMismatch[4]
			print ''

		print '=================================================================================='


	#==============================================================================================#
	def loadGenome(self, genomeFileName):

		f = open(genomeFileName)
		genomeStr = f.readline()
		genomeStr = genomeStr.strip()
		genomeStr = genomeStr[1:]
		f.close()

		return genomeStr
