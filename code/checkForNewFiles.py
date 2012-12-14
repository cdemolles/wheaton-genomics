from ftplib import FTP
from ConfigParser import SafeConfigParser
import os
import pickle
import smtplib
import re


#==============================================================================================
def main():

	# open the configuration file
	parser = SafeConfigParser()
	parser.read('config.ini')

	# holds a cached copy of the dictionary
	pickleFileName = os.path.join(parser.get('main', 'script_dir'), parser.get('downloader', 'cached_file_name'))

	# holds a list of new files to download
	newFiles = []

	# gets the cached list of genome files
	genomeFiles = getCachedList(pickleFileName)

	# directory where genome files are located
	ftpDirectory = parser.get('downloader', 'ftp_dir')

	# connect to the ftp server and go to ftpDirectory
	ftp = FTP(parser.get('downloader', 'ftp_site'))
	ftp.login()
	ftp.cwd(ftpDirectory)
	
	# get fna list, ptt list, and rnt list
	fnaFiles = ftp.nlst('*/*.fna')
	pttFiles = ftp.nlst('*/*.ptt')
	rntFiles = ftp.nlst('*/*.rnt')

	# run through all the fna files and add any ones not in the dictionary to the list of things to be downloaded
	for fnaFile in fnaFiles:
		if not genomeFiles.has_key(fnaFile):
			newFiles.append(fnaFile)

	# run through all the ptt files and add any ones not in the dictionary to the list of things to be downloaded
	for pttFile in pttFiles:
		if not genomeFiles.has_key(pttFile):
			newFiles.append(pttFile)

	# run through all the rnt files and add any ones not in the dictionary to the list of things to be downloaded
	for rntFile in rntFiles:
		if not genomeFiles.has_key(rntFile):
			newFiles.append(rntFile)

	# if there are files to be downloaded
	if newFiles != []:

		# loop through the list and download them
		for newFile in newFiles:

			path = parser.get('main', 'root_dir')

			print "Downloading file " + newFile
			downloadFile(newFile, path, ftp)
			print "Finished downloading file " + newFile

			genomeFiles[newFile] = True
			newFileName = renameFile(os.path.join(path, newFile))

			newFileNameExt = os.path.splitext(newFileName)[1]
			if newFileNameExt == parser.get('extensions', 'genome'):
				convertToSingleLine(newFileName, 1)

		# send an email indicating that there were items downloaded
		send = parser.getboolean('downloader', 'send_email')

		if send:
			sendEmail(['demolles_christopher@wheatoncollege.edu', 'mleblanc@wheatoncollege.edu'], genomeFiles)

		saveCachedList(pickleFileName, genomeFiles)

	else:

		print "No new files found. Nothing downloaded."


#==============================================================================================
def getCachedList(pickleFileName):
	
	# open and load the pickled file
	if os.path.isfile(pickleFileName):
		pickleFile = open(pickleFileName, 'rb')
		genomeFiles = pickle.load(pickleFile)
		pickleFile.close()
	else:
		genomeFiles = {}

	return genomeFiles


#==============================================================================================
def saveCachedList(pickleFileName, genomeFiles):

	# save the genome file dictionary to a file
	pickleFile = open(pickleFileName, 'wb')
	pickle.dump(genomeFiles, pickleFile)
	pickleFile.close()


#==============================================================================================
def downloadFile(genomeFile, location, ftp):

	localFile = os.path.join(location, genomeFile)

	dirName = os.path.dirname(localFile)

	if not os.path.exists(dirName):
		os.mkdir(dirName)

	f = open(localFile, 'wb')
	ftp.retrbinary('RETR ' + genomeFile, f.write)
	f.close()


#==============================================================================================
def renameFile(fileName):

	extension = os.path.splitext(fileName)[1]

	# open the file
	f = open(fileName, 'r')

	# read the first line (stripping off the newline character
	firstLine = f.readline().strip()

	# if the extension is fna
	if extension == '.fna':

		# split the string on the '|' character
		firstLine = firstLine.split('|')

		# the file name is the last entry in the list after split
		newFileName = firstLine[-1].strip()

	# for all other extensions
	else:

		# split the string on a ' - ' character
		firstLine = firstLine.split(' - ')

		# the file name is the first entry in the list after split
		newFileName = firstLine[0]
	
	# close the file
	f.close()

	# get the directory to the current file
	path = os.path.dirname(fileName)

	# remove problematic punctuation characters from the file name and add the appropriate extension back
	newFileName = re.sub(r'\W+', '_', newFileName) + extension

	# include the path with the file name
	newFileName = os.path.join(path, newFileName)

	# if the file name is not the same as the current file name
	if newFileName != fileName:

		# rename it
		os.rename(fileName, newFileName)

	return newFileName


#==============================================================================================
def convertToSingleLine(oldSequenceFile, linesToSkip):

	newSequenceFile = oldSequenceFile + '.oneline'

	inFile = open(oldSequenceFile, 'r')
	outFile = open(newSequenceFile, 'w')

	outFile.write('X')

	counter = 1
	for line in inFile:
		
		if counter > linesToSkip:
			outLine = line.strip()
			outFile.write(outLine)

		counter = counter + 1

	inFile.close()
	outFile.close()


#==============================================================================================
def sendEmail(toField, listOfNewFiles):
	
	newFiles = "\n".join(listOfNewFiles)

	sender = "Genomics <genomics@wheatoncollege.edu>"

	message = "From: " + sender + "\n"
	message = message + "To: " + ", ".join(toField) + "\n"
	message = message + "Subject: New Genome Sequences Found\n"
	message = "The following files were downloaded: \n"
	message = message + newFiles	

	smtpServer = smtplib.SMTP('outbound.wheatonma.edu', 25)
	smtpServer.sendmail(sender, toField, message)
	smtpServer.quit()


if __name__ == "__main__":
	main()
