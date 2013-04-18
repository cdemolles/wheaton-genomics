"""
-----------------------------------------------------------------------------------------------------------------
Bug.py
Programmer: Chris DeMolles

Tested with Python 2.7.3

This module defines the Bug class. The Bug class stores information about a bug. Each bug has a name, the kingdom
to which it belongs, and whether it is single cellular or multicelluar.
-----------------------------------------------------------------------------------------------------------------
"""

import os, sqlite3
from sequence import Sequence
from ConfigParser import SafeConfigParser


class Bug(object):

    #-----------------------------------------------------------------------------------------------------------------#
    def __init__(self, name):
        """
        ==============================================================================================
        SUMMARY: Initializes bug object
        PRE: Expects name to be the name of a bug in the database of bugs
        POST: Initializes a new bug object with name and genome based on that name
        ==============================================================================================
        """

        # use a parser to open the configuration file
        parser = SafeConfigParser()
        parser.read('config.ini')

        # get the data directory from the configuration file
        dataDirectory = os.path.join(parser.get('directories', 'root_dir'), parser.get('directories', 'data_subdir'))

        # get the bug directory from the configuration file
        bugDirectory  = os.path.join(dataDirectory, parser.get('directories', 'bug_subdir'))

        # get the database file name and the name of the bug
        self.databaseFileName = os.path.join(dataDirectory, parser.get('files', 'database_file_name'))
        
        # make a new connection to the database
        dbConnection = sqlite3.connect(self.databaseFileName)
        cursor = dbConnection.cursor()
        
        # query the database for the information about this bug
        bugDataQuery = "SELECT name, type, subfolder, file_prefix, kingdom, cellular, category FROM organisms WHERE name = ?"
        parameters = (name,)

        # execute the query
        cursor.execute(bugDataQuery, parameters)

        # fetch the results from this query (at most one row since name is a primary key)
        bugData = cursor.fetchone()

        # if the results are not None
        if bugData is not None:

            # retrieve the name of the file where the genome is stored
            genomeFileName = os.path.join(bugDirectory + '/' + bugData[1].encode('ascii', 'ignore') + '/' + bugData[2].encode('ascii', 'ignore') + '/' + bugData[3].encode('ascii', 'ignore') + '.fna.oneline')
            
            # initialize the following data members of bug
            self.bugName  = bugData[0].encode('ascii', 'ignore')
            self.kingdom  = bugData[4]#.encode('ascii', 'ignore')
            self.cellular = bugData[5]#.encode('ascii', 'ignore')
            self.group    = bugData[6]#.encode('ascii', 'ignore')
            self.genome   = Sequence(self.loadGenome(genomeFileName))
            
        # otherwise, no results found, so the bug name must be incorrect
        else:

            print 'Invalid bug name. Check to make sure the bug name is spelled correctly or its genome file exists in its directory.'
        
        # close the connection to the database
        dbConnection.close()
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def __str__(self):
        """
        ==============================================================================================
        SUMMARY: Overloads the string casting operator
        PRE: assigned(self.genome)
        POST: Returns self.genome.sequence
        ==============================================================================================
        """

        return self.genome.sequence
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def __len__(self):
        """
        ==============================================================================================
        SUMMARY: Returns the length of the bug's genome
        PRE: assigned(self.genome)
        POST: Returns the length of self.genome
        ==============================================================================================
        """

        return len(self.genome)
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def __getitem__(self, key):
        """
        ==============================================================================================
        SUMMARY: Overloads indexing and slicing operator
        PRE: assigned(self.genome), assigned(key)
        POST: Returns self.genome[key]
        ==============================================================================================
        """

        # indexes or slices the genome and returns it
        sequenceSlice = self.genome[key]

        return sequenceSlice
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def getSequences(self, type):
        """
        ==============================================================================================
        SUMMARY: Retrieves a dictionary of sequences of a user-defined type (RNA or DNA)
        PRE: assigned(self.bugName), assigned(self.databaseName), assigned(type)
        POST: Returns a dictionary with the starting position of the sequence as the key and the sequence as the value
        ==============================================================================================
        """

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

        # for each row in the result set
        for row in cursor.execute(sqlString, parameters):
    
            # assign the following attributes from the database values
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
            
            # create a new sequence object with the attributes above
            sequence = Sequence(sequenceStr, start, end, strand, type, name)
            sequences[start] = sequence

        # close the connection to the database
        dbConnection.close()

        return sequences
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def getUniqueRNA(self):
        """
        ==============================================================================================
        SUMMARY: Retrieves a dictionary of sequences of the unique RNA in this bug
        PRE: assigned(self.bugName), assigned(self.databaseName)
        POST: Returns a dictionary with the name of the sequence as the key and the sequence as the value
        ==============================================================================================
        """

        # initialize an empty dictionary to hold the sequences
        sequences = {}

        # connect to the database
        dbConnection = sqlite3.connect(self.databaseFileName)

        # make a new database cursor
        cursor = dbConnection.cursor()

        # store the parameters for the query
        parameters = (self.bugName, 'RNA',)

        # extract the sequence data from the database given the name of the bug and the type of sequence
        sqlString = "SELECT sequence, starting_location, ending_location, strand, type, sequence_name FROM sequences WHERE organism_name = ? AND type = ? GROUP BY sequence_name"

        # for each row in the result set
        for row in cursor.execute(sqlString, parameters):

            # assign the following attributes from the database values
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
            
            # create a new sequence object with the attributes above
            sequence = Sequence(sequenceStr, start, end, strand, type, name)
            sequences[name] = sequence

        # close the connection to the database
        dbConnection.close()

        return sequences
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def getRNA(self):
        """
        ==============================================================================================
        SUMMARY: Retrieves a dictionary of sequences of the RNA in this bug
        PRE: assigned(self.bugName), assigned(self.databaseName)
        POST: Returns a dictionary with the starting position of the sequence as the key and the sequence as the value
        ==============================================================================================
        """

        return self.getSequences('RNA')
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def getDNA(self):
        """
        ==============================================================================================
        SUMMARY: Retrieves a dictionary of sequences of integenic DNA in this bug
        PRE: assigned(self.bugName), assigned(self.databaseName)
        POST: Returns a dictionary with the starting position of the sequence as the key and the sequence as the value
        ==============================================================================================
        """

        return self.getSequences('DNA')
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def getAnnotation(self, start, end):
        """
        ==============================================================================================
        SUMMARY: Gets a list of sequences between start and end to prepare for annotation
        PRE: assigned(self.bugName), assigned(self.databaseName), assigned(start), assigned(end)
        POST: Returns a list of sequences between start and end
        ==============================================================================================
        """

        # make a blank list to hold the sequences for annotation
        annotationList = []

        # connect to the database
        dbConnection = sqlite3.connect(self.databaseFileName)

        # make a new database cursor
        cursor = dbConnection.cursor()

        # store the parameters for the query
        parameters = (self.bugName, start, end, start, end, start, end)

        # for each row in the result set
        for row in cursor.execute("SELECT sequence, starting_location, ending_location, strand, type, sequence_name FROM sequences WHERE organism_name = ? AND (((? BETWEEN starting_location AND ending_location) OR (? BETWEEN starting_location AND ending_location)) OR ((starting_location BETWEEN ? AND ?) OR (ending_location BETWEEN ? AND ?))) ORDER BY starting_location, ending_location", parameters):

            # assign the following attributes from the database values
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
            
            # create a new sequence object with the attributes above
            sequence = Sequence(sequenceStr, start, end, strand, type, name)

            # add the sequence to the list
            annotationList.append(sequence)

        # close the connection to the database
        dbConnection.close()

        return annotationList
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def shuffle(self):
        """
        ==============================================================================================
        SUMMARY: Duplicates this bug and shuffles the genome of the new bug
        PRE: assigned(self.bugName), assigned(self.databaseName), assigned(start), assigned(end)
        POST: Returns a list of sequences between start and end
        ==============================================================================================
        """

        # make a new bug object with the same information as this bug
        shuffledBug = Bug(self.bugName)

        # shuffle the new bug's genome
        shuffledBug.genome = shuffledBug.genome.shuffle()

        return shuffledBug
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def printReport(self):
        """
        ==============================================================================================
        SUMMARY: Prints a report about the integenic DNA and RNA sequences for the entire genome
        PRE: assigned(self.bugName), assigned(self.databaseName)
        POST: Prints a report about the intergenic DNA and RNA sequences to the console
        ==============================================================================================
        """

        # get the intergenic DNA and RNA sequences for this bug
        dna = self.getDNA()
        rna = self.getRNA()
        
        # print the name of the bug and the length of its genome
        print '=================================================================================='
        print 'The name of this bug is', self.name
        print 'The length of its genome is', len(self.genome)
        print ''

        # for each RNA in the genome
        for rnaName, sequence in rna.iteritems():

            # count the number of potential perfect inverted repeats of stem-length 3 and 4
            perfectIRs = sequence.countInvertedRepeats(3, 4, 0)

            # count the number of potential inverted repeats with one mismatch of stem-length 3 and 4
            IRsOneMismatch = sequence.countInvertedRepeats(3, 4, 1)

            # change +/- to direct/indirect
            if sequence.strand == '+':
                strand = 'direct'
            else:
                strand = 'indirect'

            # print the report to the console
            print 'There is a', rnaName, 'RNA on the', strand, 'strand that begins at', sequence.start, 'and ends at', sequence.end
            print 'The length of the sequence is', len(sequence)
            print 'The ratio of potential perfect IRs of length 6 to the length of the sequence is', perfectIRs[3]
            print 'The ratio of potential perfect IRs of length 8 to the length of the sequence is', perfectIRs[4]
            print 'The ratio of potential IRs of length 6 with 1 mismatch to the length of the sequence is', IRsOneMismatch[3]
            print 'The ratio of potential IRs of length 8 with 1 mismatch to the length of the sequence is', IRsOneMismatch[4]
            print ''

        # for each intergenic DNA in the genome
        for dnaStart, sequence in dna.iteritems():

            # count the number of potential perfect inverted repeats of stem-length 3 and 4
            perfectIRs = sequence.countInvertedRepeats(3, 4, 0)

            # count the number of potential inverted repeats with one mismatch of stem-length 3 and 4
            IRsOneMismatch = sequence.countInvertedRepeats(3, 4, 1)

            # change +/- to direct/indirect
            if sequence.strand == '+':
                strand = 'direct'
            else:
                strand = 'indirect'

            # print the report to the console
            print 'There is a DNA sequence on the', strand, 'strand that begins at', sequence.start, 'and ends at', sequence.end
            print 'The length of the sequence is', len(sequence)
            print 'The ratio of potential perfect IRs of length 6 to the length of the sequence is', perfectIRs[3]
            print 'The ratio of potential perfect IRs of length 8 to the length of the sequence is', perfectIRs[4]
            print 'The ratio of potential IRs of length 6 with 1 mismatch to the length of the sequence is', IRsOneMismatch[3]
            print 'The ratio of potential IRs of length 8 with 1 mismatch to the length of the sequence is', IRsOneMismatch[4]
            print ''

        print '=================================================================================='
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def loadGenome(self, genomeFileName):
        """
        ==============================================================================================
        SUMMARY: Loads the bug's genome from genomeFileName
        PRE: assigned(genomeFileName)
        POST: Stores the bug's genome in a string
        ==============================================================================================
        """

        # open genomeFileName
        f = open(genomeFileName)

        # read the contents of the file and store it in genomeStr
        genomeStr = f.readline()

        # strip any whitespace from the string
        genomeStr = genomeStr.strip()

        # remove the first character from the string (there is a bogus 'X' character at the beginning of the string)
        genomeStr = genomeStr[1:]

        # close the file
        f.close()

        return genomeStr
    #-----------------------------------------------------------------------------------------------------------------#
