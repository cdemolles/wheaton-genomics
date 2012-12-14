Step 1: Edit config.ini
Step 2: Run checkForNewFiles.py
	Outcome: We now have all the bugs that have legitimate ptt, fna, and rnt files (all other bugs are skipped)
Step 3: Assuming bug.py and sequence.py is available, run motifCounter.py
	Outcome: Produces SVM_DNA.txt and SVM_RNA.txt, both in preparation for SVM
Step 4: Assuming SVM_DNA.txt and SVM_RNA.txt are available, run dna_rna.py
	Outcome: Runs scilearnkit with combined data
