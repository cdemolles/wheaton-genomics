#Wheaton College Genomics Research Group

This software repository contains source code for our research project to better understand the underlying structures of RNA.

##Documentation

###Organization
The ``code`` directory contains the source code for this project. The ``data`` directory contains the bug data from the HMP website and a SQLite database that contains the the RNA and intergenic DNA sequences for each bug.

###Setting up the repository on your machine
1. Clone the repository to your machine by doing a ``git clone``.
2. Descend into the ``code`` directory.
3. Open ``config.ini`` for editing.
4. Edit the ``root_dir`` entry under ``[directories]`` to match the location of where the repository is stored on your machine.

###Running the code
After performing the steps above, you may now run ``count.py`` in the ``code`` directory. ``count.py`` expects that the SQLite database ``genomics.sqlite`` is located in the ``data`` directory in order to run properly.

``count.py`` exports a tab-separated value file containing the RNA and intergenic DNA sequences for each bug. The counts of each of the 256 4-mers, the number of potential perfect IRs of stem-length 3,4,5, and the number of potential IRs of stem-length 3,4,5 with one mismatch are saved to the file. The tab-separated value file is saved as ``IR_counts_345.tsv`` in the ``data`` directory.