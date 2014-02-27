#!/bin/bash

#####################Making Directories###############################
echo Making MFS and Output directories if they have not already been created...
mkdir -p MFS
mkdir -p MFS/Output/
#copying below file to MFS directory. This is a hack to make sure error message 
#does not occur when running the conditional statements below. This script can be
#run if this file is not copied now, but you will get an error message until the
#MFS directory is not empty.
cp 93279115_MFS.gb MFS/93279115_MFS.gb

#####################Calling getGI.py script###############################
echo Downloading files...
python getGI.py

####################File size checking###############################
#This script checks to make sure all files are < 2kb. If so, the file
#is deleted, re-downloaded and then re-run. Script works recrussively.
Accessioncount=`awk 'END {print NR}' accessionList.txt` #or RequiredNumber = $(wc -l accessionList.txt)
Filecount=$(find MFS/*.gb -type f -print| wc -l)
if [[ $Filecount -eq $Accessioncount ]]; then
	echo 
	echo Checking for corrupted files...
	echo 
	python getSize.py
fi

#####################Copying 93279115_MFS.gb to MFS Directory###############################
#The original file has errors in it. This file must be copied to MFS/ otherwise NextractGb.py will not function!
#Don't even think about not copying this file over to the directory.
Filecount=$(find MFS/*.gb -type f -print| wc -l)
if [[ $Filecount -eq $Accessioncount ]]; then
echo Copying file 93279115_MFS.gb...
	echo 
	cp 93279115_MFS.gb MFS/93279115_MFS.gb
fi

################################Verifying number of files####################################
Filecount=$(find MFS/*.gb -type f -print| wc -l)

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
