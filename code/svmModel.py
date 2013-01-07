import numpy as np
import pylab as pl
from sklearn import datasets, svm
from sklearn.datasets import load_svmlight_file

def main():
	
	labelsRNA = np.genfromtxt("/home/chris/SVM_RNA.txt", dtype=str, skip_header=1, usecols=(0))
	dataRNA = np.genfromtxt("/home/chris/SVM_RNA.txt", dtype=float, skip_header=1, usecols=tuple(range(257,261)), converters={260:classificationToInt})

	labelsDNA = np.genfromtxt("/home/chris/SVM_DNA.txt", dtype=str, skip_header=1, usecols=(0))
	dataDNA = np.genfromtxt("/home/chris/SVM_DNA.txt", dtype=float, skip_header=1, usecols=tuple(range(257,261)), converters={260:classificationToInt})
	
	# extract the last column from every row (i.e. the classification for the sequence -- 0 == notRNA, 1 == RNA)
	targetRNA = dataRNA[:,-1]
	targetDNA = dataDNA[:,-1]

	# leave out the classification column so the data can be used as a dataset
	datasetRNA = dataRNA[:, :-1]
	datasetDNA = dataDNA[:, :-1]

	# get the size of the sample, i.e. the length of dataset
	sampleSizeRNA = len(datasetRNA)
	sampleSizeDNA = len(datasetDNA)
	
	# define the percentage of the dataset to be used for training
	percentage = .9
	testing = 250

	# see random generator
	np.random.seed()

	# create a permutation of sampleSize numbers (to be used to shuffle the dataset and the target)
	orderRNA = np.random.permutation(sampleSizeRNA)
	orderDNA = np.random.permutation(sampleSizeDNA)
	
	# shuffle the dataset, target, and the labels
	datasetRNA = datasetRNA[orderRNA]
	targetRNA  = targetRNA[orderRNA].astype(np.int)
	labelsRNA  = labelsRNA[orderRNA]

	datasetDNA = datasetDNA[orderDNA]
        targetDNA  = targetDNA[orderDNA].astype(np.int)
        labelsDNA  = labelsDNA[orderDNA]

	# training set (90% of the sample)
	# take the first percentage * sampleSize rows along with all the columns from the dataset
	datasetTrainRNA = datasetRNA[:testing]
	targetTrainRNA  = targetRNA[:testing]
	labelsTrainRNA  = labelsRNA[:testing]

	datasetTrainDNA = datasetDNA[:testing]
        targetTrainDNA  = targetDNA[:testing]
        labelsTrainDNA  = labelsDNA[:testing]

	# testing set (10% of the sample)
	# take the last percentage * sampleSize rows along with all the columns from the dataset
	datasetTestRNA = datasetRNA[testing:]
	targetTestRNA  = targetRNA[testing:]
	labelsTestRNA  = labelsRNA[testing:]

	datasetTestDNA = datasetDNA[testing:]
        targetTestDNA  = targetDNA[testing:]
        labelsTestDNA  = labelsDNA[testing:]

	datasetTrain = np.concatenate((datasetTrainRNA, datasetTrainDNA))
	targetTrain  = np.concatenate((targetTrainRNA, targetTrainDNA))
	labelsTrain  = np.concatenate((labelsTrainRNA, labelsTrainDNA))

	datasetTest  = np.concatenate((datasetTestRNA, datasetTestDNA))
	targetTest   = np.concatenate((targetTestRNA, targetTestDNA))	
	labelsTest   = np.concatenate((labelsTestRNA, labelsTestDNA))

	
	
	# build the SVM model	
	clf = svm.SVC(kernel='rbf', gamma=10)
	
	# fit the data
	clf.fit(datasetTrain, targetTrain)

	# predict the target values for datasetTest using the model
	results = clf.predict(datasetTest).astype(np.int)

	# count how many were labeled correctly
	countsCorrect   = {'RNA':0, 'notRNA':0}

	# count how many RNA were labeled as notRNA and how many notRNA were labeled as RNA
	countsIncorrect = {'RNA':0, 'notRNA':0}

	print 'Label\tPredicted Classification\tActual Classification'	
	for i in range(len(results)):

		print labelsTest[i] + '\t' + str(results[i]) + '\t' + str(targetTest[i])

		if results[i] == targetTest[i]:
			countsCorrect[intToClassification(targetTest[i])] = countsCorrect[intToClassification(targetTest[i])] + 1
		else:
			countsIncorrect[intToClassification(targetTest[i])] = countsIncorrect[intToClassification(targetTest[i])] + 1

	print 'RNA labeled correctly:    ', countsCorrect['RNA']
	print 'RNA mislabeled as notRNA: ', countsIncorrect['RNA']
	print 'notRNA labeled correctly: ', countsCorrect['notRNA']
	print 'notRNA mislabeled as RNA: ', countsIncorrect['notRNA']

def classificationToInt(x):

	if x == 'notRNA':
		return 0.0
	else:
		return 1.0

if __name__ == '__main__':
	main()
