from bug import *
from sequence import *
import sqlite3
import os
import re
import dnaFunctions

def main():

    parser = SafeConfigParser()
    parser.read('config.ini')

    # get the root directory from the configuration file
    rootDirectory = parser.get('directories', 'root_dir')

    # get the data directory from the configuration file
    dataDirectory = os.path.join(rootDirectory, parser.get('directories', 'data_subdir'))

    # get the database file name
    databaseFileName = os.path.join(dataDirectory, parser.get('files', 'database_file_name'))

    dbConnection = sqlite3.connect(databaseFileName)
    cursor = dbConnection.cursor()

    report = open('../data/report.tsv', 'w')
    outputHeading(report)

    bugDataQuery = "SELECT name, kingdom, category FROM organisms WHERE kingdom <> 'Human_Microbiom' AND kingdom <> 'Eukaryotic'"

    # for each tuple of bug data in the result set
    for bugData in cursor.execute(bugDataQuery):

        # the bug name is the first column in the result set
        bugName  = bugData[0].encode('ascii', 'ignore')
        kingdom  = bugData[1].encode('ascii', 'ignore')
        category = bugData[2]

        if category is not None:
            category = category.encode('ascii', 'ignore')
        else:
            category = ''

        print "Counting", bugName, "..."

        # make a new bug object given the bug name
        bug = Bug(bugName)

        # gets a dictionary of RNA sequences (with the key as the name of the RNA and the value as a sequence object)
        bugRNA = bug.getRNA()

        # gets a dictionary of DNA sequences (with the key as the starting location of the DNA and the value as a sequence object)
        bugDNA = bug.getDNA()

        outputData(report, bugRNA, bugName, kingdom, category, 'RNA')
        outputData(report, bugDNA, bugName, kingdom, category, 'DNA')

        print "Done counting", bugName, "..."

    report.close()


def outputHeading(outputFile):

    outputFile.write('Sequence_Name')
    outputFile.write('\t')

    outputFile.write('Bug_Name')
    outputFile.write('\t')

    outputFile.write('Kingdom')
    outputFile.write('\t')

    outputFile.write('Category')
    outputFile.write('\t')

    outputFile.write('Type')
    outputFile.write('\t')

    outputFile.write('Sequence_Length')
    outputFile.write('\t')

    outputFile.write('IRs_Stem_Len_3_Perfect')
    outputFile.write('\t')

    outputFile.write('IRs_Stem_Len_4_Perfect')
    outputFile.write('\t')

    outputFile.write('IRs_Stem_Len_5_Perfect')
    outputFile.write('\t')

    outputFile.write('Total_IRs_Stem_Len_3')
    outputFile.write('\t')

    outputFile.write('Total_IRs_Stem_Len_4')
    outputFile.write('\t')

    outputFile.write('Total_IRs_Stem_Len_5')   
    outputFile.write('\n')


def outputData(outputFile, data, bugName, kingdom, category, type):

    for sequenceName, sequence in data.iteritems():

        #motifs                  = sequence.countMotifs(4, 4)
        invertedRepeats         = sequence.countInvertedRepeats(3, 5, 0, False)
        totalInvertedRepeats    = sequence.countTotalInvertedRepeats(3, 5, 0)
        #invertedRepeatsMismatch = sequence.countInvertedRepeats(4, 4, 1)

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

        outputName = str(bugName) + '_' + str(type) + '_' + str(outputName) + '_' + str(sequence.end) + '_' + str(strandName)

        length = abs(sequence.end - sequence.start) + 1

        outputFile.write(outputName)
        outputFile.write('\t')

        outputFile.write(bugName)
        outputFile.write('\t')

        outputFile.write(kingdom)
        outputFile.write('\t')

        outputFile.write(category)
        outputFile.write('\t')

        if type == 'DNA':
            outputFile.write('Intergenic')
        else:
            outputFile.write(type)

        outputFile.write('\t')

        outputFile.write(str(length))
        outputFile.write('\t')

        outputFile.write(str(invertedRepeats[3]))
        outputFile.write('\t')

        outputFile.write(str(invertedRepeats[4]))
        outputFile.write('\t')

        outputFile.write(str(invertedRepeats[5]))
        outputFile.write('\t')

        outputFile.write(str(totalInvertedRepeats[3]))
        outputFile.write('\t')

        outputFile.write(str(totalInvertedRepeats[4]))
        outputFile.write('\t')

        outputFile.write(str(totalInvertedRepeats[5]))
        outputFile.write('\n')
        

if __name__ == '__main__':
    main()
