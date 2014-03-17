#!/bin/bash

aligned=$1
accession=$2
db=$3

echo db

#####################Making Directories###############################
echo Making $db and Output directories if they have not already been created...
mkdir -p $db
mkdir -p $db/Output/

#copying below file to $db directory. This is a hack to make sure error message 
#does not occur when running the conditional statements below. This script can be
#run if this file is not copied now, but you will get an error message until the
#$db directory is not empty.
cp 93279115_$db.gb $db/93279115_$db.gb


#####################Calling getGI.py script###############################
echo Downloading files...

python getGI.py $aligned $accession $db

####################File size checking###############################
#This script checks to make sure all files are < 2kb. If so, the file
#is deleted, re-downloaded and then re-run. Script works recrussively.
Accessioncount=`wc -l accessionList.txt`
Filecount=`find $db/*.gb -type f | wc -l`
if [[ $Filecount -eq $Accessioncount ]]; then
	echo 
	echo Checking for corrupted files...
	echo 
	python getSize.py
fi

#TODO wtf is this 

# Copying 93279115_$db.gb to $db Directory
# The original file has errors in it. This file must be copied to $db/ otherwise NextractGb.py will not function!
# Don't even think about not copying this file over to the directory.
Filecount=$(find $db/*.gb -type f -print| wc -l)
if [[ $Filecount -eq $Accessioncount ]]; then
echo Copying file 93279115_$db.gb...
	echo 
	cp 93279115_$db.gb $db/93279115_$db.gb
fi

# Verifying number of files
Filecount=$(find $db/*.gb -type f | wc -l)

if [[ $Filecount -eq $Accessioncount ]]; then
	echo Extracting data...
	python NextractGB.py
	echo Looks like you are ready to look at some data!
	
else
	echo 
	echo There are files that need to be downloaded or files that are corrupted and must be re-downloaded. 
	echo If you are downloading more than 100 files/day, files can only be downloaded between 9pm-5am EST on weekdays. 
	echo If it is not the weekend and you are outside of those hours, try running this again at an appropriate hour.
fi
