"""
-----------------------------------------------------------------------------------------------------------------
Sequence.py
Programmer: Chris DeMolles

Tested with Python 2.7.3

This module defines the Sequence class. The Sequence class stores information about a sequence. Each sequence has
a starting location, ending location, the strand it's on, the type of sequence, and the name of the sequence.

Methods include generating the motifs for the sequence and counting the number of inverted repeats for the sequence.
-----------------------------------------------------------------------------------------------------------------
"""

import dnaFunctions as DNA
import random, re


class Sequence(object):

    #-----------------------------------------------------------------------------------------------------------------#
    def __init__(self, sequence, start=None, end=None, strand=None, type=None, name=None):
        """
        ==============================================================================================
        SUMMARY: Initializes sequence object
        PRE: assigned(sequence), optional(start), optional(end), optional(strand), optional(type), optional(name)
        POST: Initializes a new sequence object with the provided parameters
        ==============================================================================================
        """

        # assign the class members the following parameter data
        self.sequence = sequence.strip()
        self.start    = start
        self.end      = end
        self.strand   = strand
        self.type     = type
        self.name     = name
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def __str__(self):
        """
        ==============================================================================================
        SUMMARY: Overloads the string casting operator
        PRE: assigned(self.sequence)
        POST: Returns self.sequence
        ==============================================================================================
        """

        return self.sequence
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def __len__(self):
        """
        ==============================================================================================
        SUMMARY: Returns the length of the sequence
        PRE: assigned(self.sequence)
        POST: Returns the length of self.sequence
        ==============================================================================================
        """

        return len(self.sequence)
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def __getitem__(self, key):
        """
        ==============================================================================================
        SUMMARY: Overloads indexing and slicing operator
        PRE: assigned(self.sequence), assigned(key)
        POST: Returns a new sequence initialized to self.sequence[key]
        ==============================================================================================
        """

        # add bogus character to beginning of sequence for 1-based indexing
        section = 'X' + self.sequence

        return Sequence(section[key])
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def reverseComplement(self):
        """
        ==============================================================================================
        SUMMARY: Returns the reverse complement of the sequence
        PRE: assigned(self.sequence)
        POST: Returns the reverse complement of self.sequence
        ==============================================================================================
        """

        # find the reverse complement of self.sequence
        reverseComplement = DNA.reverseComplement(self.sequence)

        # convert the reverse complement to a sequence object
        reverseComplement = Sequence(reverseComplement)

        return reverseComplement
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def motifs(self, min, max):
        """
        ==============================================================================================
        SUMMARY: Finds all the motifs of length min to length max (inclusive)
        PRE: assigned(self.sequence), assigned(min), assigned(max)
        POST: Returns a dictionary with the key as the length of the motif and the value as each of
              the motifs of that length
        ==============================================================================================
        """

        # tack on a dummy character to the sequence so the actual first character in the sequence
        # starts at index 1
        section = 'X' + self.sequence

        # get the length of the sequence
        strandLength = len(section)

        # start at the first character in the sequence
        start  = 1

        # initialize a blank dictionary to store all motifs of length min up to length max
        motifs = {}
    
        # while haven't traversed the entire genome
        while start < len(section):
        
            # for i from min to max (inclusive)
            for i in xrange(min, max + 1): 

                # if the start is less than the length of the strand - i
                if start <= strandLength - i:

                    # cut from start to start + i to form the motif
                    motif = section[start:start + i]

                    # if the motif is already in the dictionary
                    if motif in motifs:
                        # append the starting location to the dictionary of motifs
                        motifs[motif].append(start)
                    # otherwise
                    else:
                        # initialize the list with start
                        motifs[motif] = [start]
        
            # increment start
            start = start + 1
    
        return motifs
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def getInvertedRepeats(self, min, max, mismatches, type='pairings'):
        """
        ==============================================================================================
        SUMMARY: Finds potential inverted repeats of stem-length min to max (inclusive) with mismatches
        PRE: assigned(self.sequence), assigned(min), assigned(max), assigned(mismatches) (0 or 1)
        POST: Returns a dictionary with the key as the stem-length of the IR and the value as the
              number of IRs of that stem-length
        ==============================================================================================
        """

        # initialize a blank dictionary to hold the inverted repeats
        irs = {}

        # get the motifs from min to max
        motifs = self.motifs(min, max)

        # loop over each motif in the list of motifs
        for motif in motifs:

            # make a set to store the potential inverted repeats
            # using a set to avoid adding the same element multiple times
            invertedRepeats = set()

            # take the reverse complement of the motif
            reverseComplementMismatches = DNA.reverseComplementMismatches(motif, mismatches)

            # for each reverse complement with one or more mismatched base pairs
            for reverseComplement in reverseComplementMismatches:

                # if the reverse complement is also in the dictionary, then we know we have a possible IR
                if reverseComplement in motifs:

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

                                # depending on what we want to report, add different things to the inverted repeats set
                                if type == 'pairings':

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
                                        invertedRepeats.add(identifier)

                                elif type == 'nucleotides':

                                    # get the locations of the nucleotides that make up this invertedRepeat
                                    motifSpace = set(range(motifLocation, motifLocation + len(motif)))
                                    reverseComplementSpace = set(range(reverseComplementLocation, reverseComplementLocation + len(reverseComplement)))

                                    # add them to the set
                                    invertedRepeats.update(motifSpace)
                                    invertedRepeats.update(reverseComplementSpace)                                    

            # get the length of the motif
            key = len(motif)

            # if we've already begun counting IRs of length 2 * key
            if key in irs:
                # union the exisiting set with the generated set of inverted repeats
                irs[key] = irs[key].union(invertedRepeats)
            # otherwise
            else:
                # initialize and empty set
                irs[key] = set()

        return irs
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def countMotifs(self, min, max, relativePercentages=False):
        """
        ==============================================================================================
        SUMMARY: Counts the number of motifs of length min to length max (inclusive)
        PRE: assigned(self.sequence), assigned(min), assigned(max), optional(relativePercentages)
        POST: Returns a dictionary with the key as the length of the motif and the value as the number
               of motifs of that length
        ==============================================================================================
        """

        # initialize a blank dictionary to hold motif counts
        motifCounts = {}

        # get the motifs from min to max
        motifs = self.motifs(min, max)

        # for each motif in motifs
        for motif in motifs:

            # the number of motifs "motif" is the length of the list stored in motifs[motif]
            motifCounts[motif] = len(motifs[motif])

            # if we want relative percentages
            if relativePercentages == True:
                # divide the count by the length of the sequence
                motifCounts[motif] = float(motifCounts[motif]) / len(self.sequence)

        return motifCounts
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def countInvertedRepeats(self, min, max, mismatches, relativePercentages=False, type='pairings'):
        """
        ==============================================================================================
        SUMMARY: Counts the number of potential inverted repeats from stem-length min to length max (inclusive)
        PRE: assigned(self.sequence), assigned(min), assigned(max), assigned(mismatches), optional(relativePercentages), optional(type)
        POST: Returns a dictionary with the key as the length of the stem and the value as the number
               of IRs of that stem-length
        ==============================================================================================
        """

        # get the inverted repeats from min to max with mismatches of type
        irIdentifiers = self.getInvertedRepeats(min, max, mismatches, type)

        # initialize a blank dictionary to hold the counts
        countsOfIRs = {}

        # for each key in the inverted repeats
        for key in irIdentifiers:

            # the number of IRs of that key is the length of the list stored in irIdentifiers[key]
            countsOfIRs[key] = len(irIdentifiers[key])

            # if we want relative percentages
            if relativePercentages == True:
                # divide the count by the length of the sequence
                countsOfIRs[key] = float(countsOfIRs[key]) / len(self.sequence)

        # add 0's to the dictionary for any inverted repeat that does not appear
        for i in xrange(min, max+1):
            if i not in countsOfIRs:
                countsOfIRs[i] = 0

        return countsOfIRs
    #-----------------------------------------------------------------------------------------------------------------#


    #-----------------------------------------------------------------------------------------------------------------#
    def shuffle(self):
        """
        ==============================================================================================
        SUMMARY: Shuffles the sequence
        PRE: assigned(self.sequence)
        POST: Returns a new sequence with the shuffled version of this sequence
        ==============================================================================================
        """

        # convert the sequence to a list of characters
        listOfChars = list(self.sequence)

        # seed random
        random.seed()

        # for i from len(listOfCharacters) - 1 down to 0
        for i in range(len(listOfChars) - 1, 0, -1):

            # generate a random integer from 0 to i
            j = random.randint(0, i)

            # swap the ith and jth characters in listOfChars
            temp = listOfChars[j]
            listOfChars[j] = listOfChars[i]
            listOfChars[i] = temp

        # convert the list back to a string
        shuffledSequence = ''.join(listOfChars)

        return Sequence(shuffledSequence)
    #-----------------------------------------------------------------------------------------------------------------#
