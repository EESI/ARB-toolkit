'''
Created: E. Reichenberger (err29@drexl.edu)
Date: 9.20.2012

Purpose: Interruptions during file downloading process occur and files are not
complete. Script checks file size, deletes it when less than 2bytes, and
re-downloads file. Works recursively.

Steps:
1. Get list of all files in db directory. 
2. Check size of all files in list
3. If file is smaller than 2b, add name to new list and delete file
	a. remove db and '_db.gb' from file name
4. Re-download all items in new list.
5. Re-run script until there are no files to download.
	a. Using GI number in call
	b. IMPORTANT!!!!! If more than 100 files need to be downloaded, this can only
		 occur between 9pm-5am EST only (M-F).

IMPORTANT!!!!! Times are listed for EST, if you are in another time zone, make
the proper hour changes in this script

As of Feb 15, 2012, the default download file type is in XML format. Add
retmode='text' to keep in in gb format

e.g. handle=Entrez.efetch(db='nucleotide',id=gi[i], rettype='gb', retmode='text')
'''

import os
import glob
import sys
import re #regular expressions
from Bio import Entrez
import datetime
import time
from urllib2 import HTTPError, URLError
import httplib

fileList = []
sizeList = []

db = sys.argv[1]

#########################################################################
# Get list of files already downloaded from NCBI
#########################################################################
path = db
for infile in glob.glob( os.path.join(path, '*.gb') ): #file extension type
	fileList.append(infile)

for files in fileList:
	filesize = os.path.getsize(files)
	newFile = files.replace(db + '/' '').replace('_' + db + '.gb', '')
	if filesize < 2.0:
		sizeList.append(newFile)
		os.remove(files) #technically, this is not necessary. New downloaded files will overwrite existing ones.
sizeList.sort()

dupSizeList = []
for dups in sizeList:
	newName = db + '/' + dups + '_' + db + '.gb'
	dupSizeList.append(newName)
	
while len(sizeList) > 0:
	#########################################################################
	# Remove items in fileList from accession list
	#########################################################################

	sizeList.sort(lambda x,y: cmp(len(y), len(x)))

	print 'You have ' + str(len(sizeList)) + ' file(s) to download.'

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
	Entrez.email = 'err29@drexel.edu'
	start_day = datetime.date.today().weekday() # 0 is Monday, 6 is Sunday
	start_time = datetime.datetime.now().time()

	if ((start_day < 5 and start_time > datetime.time(hour=21)) or (start_day < 5 and start_time < datetime.time(hour=5)) or start_day > 5 or len(sizeList) <= 100 ):
		print 'Calling Entrez...'
		for gi in sizeList:
			if ( (datetime.date.today().weekday() < 5 and datetime.datetime.now().time() > datetime.time(hour=21)) or (datetime.date.today().weekday() < 5 and datetime.datetime.now().time() < datetime.time(hour=5)) or (datetime.date.today().weekday() == start_day + 1 and datetime.datetime.now().time() < datetime.time(hour=5)) or (datetime.date.today().weekday() > 5) or len(sizeList) <= 100 ):
				while True:
					print 'Downloading ' + gi					
					try:
						handle=Entrez.efetch(db='nucleotide',id=gi, rettype='gb', retmode='text') 
						FILENAME = db + '/' + gi + '_' + db + '.gb'
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
				break
	else:
		print 'You cannot run the script right now. Try again at an appropriate time.'
		break

	for size in sizeList:
		newName = db + '/' + size + '_' + db + '.gb'
		if os.path.getsize(newName) >= 2:
			sizeList.remove(size)
