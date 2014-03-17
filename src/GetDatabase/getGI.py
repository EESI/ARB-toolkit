#!/usr/bin/env python
'''
Created: E. Reichenberger (err29@drexl.edu)
Date: 7.30.2012

Purpose: Given a fasta file, extract GI information from file, place into a
list and use this information to download benbank records from Entrez. 

Notes: 
1. The file is a fasta file containing one or more sequences. 
2. Each sequence preceded by a line that starts with '>' and contains the GI.
3. The files downloaded from Entrez will be the named as follows: GInumber_MFS.gb.
4. The files will be saved to a directory (MFS) which must be created before
	 running this script (i.e. mkdir MFS).
5. In the event that the downloading process is interrupted, this script will
	 verify that a file has not already been downloaded. Script checks the MFS
   directory for *.gb files and places the file name into a list. The file name
   list is compared to the GI list; Existing filenames are removed from the GI
   list to avoid redundant downloading. 
6.'GI' and 'accession' used interchangeably.  
7. This script also creates a file (accession.txt). The number of lines in this
	 file are used as a reference number in the bash script run_script.sh

Steps:
1. Open fasta file (e.g. MFS_Align.fasta) 
2. Go through file line by line; if the line starts with '>', grab only the gi number
3. Add GI number to list and sort it alphabetically 
4. Get list of files in MFS/ directory (and remove 'MFS' and '_MFS.gb' from file name)
5. Download genbank files associated with GI number from Entrez
	a. Using GI number in call
	b. IMPORTANT!!!!! If more than 100 files need to be downloaded, this can only occur between 9pm-5am EST only (M-F).
6. Re-run script as needed if all files have not been downloaded.

Example Line:
>gi|260752245

As of Feb 15, 2012, the default download file type is in XML format. Add retmode='text' to keep in in gb format

e.g. handle=Entrez.efetch(db='nucleotide',id=gi[i], rettype='gb', retmode='text') 

IMPORTANT!!!!! Times are listed for EST, if you are in another time zone, make the proper hour changes in this script 
 


 Call script this way:

 getGI.py MFS_Align.fasta MFS_accession_list.txt MFS

'''

import os
import os.path
import glob
import sys
import re #regular expressions
from Bio import Entrez
import datetime
import time
from urllib2 import HTTPError, URLError
import httplib

# GI/accession number (from fasta file)  needed when calling Entrez
accession = [] 
fileList = []
# getAccessionList.py
# parse out the accession list from the algined fasta. Ex:
# get_accession_list.py aligned_fasta_fn accession_list_fn

# arguments
aligned_fasta_fn = sys.argv[1]
accession_list_fn = sys.argv[2]
db = sys.argv[3]

# Open file containing GI values
aligned_fasta_fh = open(aligned_fasta_fn, 'r')

# Open output file of gi numbers
if os.path.exists(accession_list_fn):
	os.remove(accession_list_fn)
accession_list_fh = open(accession_list_fn, 'w')

for line in aligned_fasta_fh:
	if line.startswith('>'):
		giLine = line.rsplit('gi|')
		giL = giLine[1].rsplit('|')
		gi = giL[0].rsplit(' ')
		accession.append(gi[0])
		accession_list_fh.write(gi[0] + '.gb\n')

accession_list_fh.close()
aligned_fasta_fh.close()


#########################################################################
# Inconsistencies exist in format of first line.
# >gi|VALUE|name|:num1 - num2 description
# >gi|VALUE:num1 - num2 description
# SEPARATION ORDER IS IMPORTANT
#########################################################################

#order list alphabetically
accession.sort()

# make sure our folder exists
if not os.path.exists(db):
    os.makedirs(db)
else:
	if not os.path.isdir(db):
		print db + "/ is not a directory"  

# Get list of files already downloaded from NCBI

# get all .gb files in our db folder
fileList = [ fn for fn in glob.glob( os.path.join(db, '*.gb')) ]

# remove basename and extension
fileList = [ os.path.basename(fn).rsplit('.', 1)[0] for fn in fileList ]

# sort our list
fileList.sort()

# Remove items in fileList from accession list
for fn in fileList:
	if fn in accession:
		accession.remove(fn)

accession.sort(lambda x,y: cmp(len(y), len(x)))

print 'You have ' + str(len(accession)) + ' file(s) to download.'

###############################################################################
# Call Entrez to download files
# If downloading more than 100 files...
# Run this script only between 9pm-5am Monday - Friday EST
# Send E-utilities requests to http://eutils.ncbi.nlm.nih.gov
# Make no more than 3 requests every 1 second.
# Use URL parameter email & tool for distributed software
# NCBI's Disclaimer and Copyright notice must be evident to users of your service. 
#
# Use this script at your own risk. 
# Neither the script author or EESI Laboratory is not responsible for consequences arising from inproper usage 
###############################################################################
start_day = datetime.date.today().weekday() # 0 is Monday, 6 is Sunday
start_time = datetime.datetime.now().time()

if ((start_day < 5 and start_time > datetime.time(hour=21)) or (start_day < 5 and start_time < datetime.time(hour=5)) or start_day > 5 or len(accession) <= 100 ):
	#print 'Calling Entrez...'
	for gi in accession:
		if ( (datetime.date.today().weekday() < 5 and datetime.datetime.now().time() > datetime.time(hour=21)) or (datetime.date.today().weekday() < 5 and datetime.datetime.now().time() < datetime.time(hour=5)) or (datetime.date.today().weekday() == start_day + 1 and datetime.datetime.now().time() < datetime.time(hour=5)) or (datetime.date.today().weekday() > 5) or len(accession) <= 100 ):
			while True:
				print 'Downloading ' + gi
				try:
					handle=Entrez.efetch(db='nucleotide',id=gi, rettype='gb', retmode='text') 
					FILENAME = db +'/' + gi + '_' + db + '.gb'
					local_file=open(FILENAME,'w')
					local_file.write(handle.read())
					handle.close()
					local_file.close()

					break
				except HTTPError, e:
					print 'Error code: ', e.code
					pass
				except URLError, e:
					print 'Reason: ', e.reason
					pass
				except httplib.HTTPException, e:
					print 'Bad Status Line'
					raise
		else:
			print 'You cannot run the script right now. Try again at an appropriate time.'
else:
	print 'You cannot run the script right now. Try again at an appropriate time.'
