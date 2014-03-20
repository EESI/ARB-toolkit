#!/usr/bin/env python
# Program: build_ift_from_metalabels.py
# Rev: 1.0
# Date: 1/25/12
# Author: S. Essinger
# Contributor: Calvin Morrison
#
# Description: This script generates the ARB .ift filter for your custom
# database. You must provide an input file containing the names of each
# meta-data field residing in your database, one name per line
# (see example: metaLables.txt). The output .ift file should then be placed
# in the ARB directory under: /arb/lib/import/                                
#
# 
# Example:
#
# build_ift_from_metalabels.py -i meta_label.txt -o output_filter.ift

import sys
import argparse

def main():
  inFileName = ""
  noLines = 0
  metaLabels = ""
  outFileName = ""

  if len(sys.argv) != 1:
    parser = argparse.ArgumentParser(description="arb IFT generator for custom databases\n")
    parser.add_argument("-i", "--input", help="input metalabel file", required=True)
    parser.add_argument("-o", "--output", help="output ift filter file", required=True)
    args = parser.parse_args()

    inFileName = args.input
    outFileName = args.output
  else:
    inFileName = raw_input('Enter name of file containing meta-labels: ')
    outFileName= raw_input('Enter desired name of ARB input filter: ')

  metaLabels = open(inFileName,'r') # Open input file for meta-label extraction
  outFile = open(outFileName,'w')

  # Include required ARB code to initialize filter.
  outFile.write('AUTODETECT      ">*"'+'\n\n')
  outFile.write('#Global settings:'+'\n')
  outFile.write('KEYWIDTH        1' + '\n\n')
  outFile.write('BEGIN           ">??*"'+'\n\n')
  outFile.write('MATCH		">*"'+'\n')
  outFile.write(r'	    SRT "* *=*1:*\t*=*1"'+'\n')
  outFile.write('	    WRITE "name"'+'\n\n')
  outFile.write('MATCH		">*"'+'\n')
  outFile.write(r'	    SRT "* *=*1:*\t*=*1"'+'\n')
  outFile.write('	    WRITE "acc"'+'\n\n')

# Create ARB code to identify, extract and assign meta-data from custom database and
# to corresponding meta-label.
  count = 1
  SRC = []
  syntax = r'"*\t*'
  syntax1 = '"* *'
  i = 0;

  lastPart = r':*\t*=*1"'
  for label in metaLabels:
    if i == 0:
      SRC.append(syntax1 + '=*' + str(count) + lastPart)
    else:
      SRC.append(syntax + '=*' + str(count) + lastPart)

    outFile.write('MATCH' + '\t' + '\t' + r'">*"' + '\n')
    outFile.write('\t' + '\t' + 'SRT ' + SRC[i] + '\n')
    outFile.write('\t' + '\t' + 'WRITE "' + str(label.strip()) + '"\n\n')

    if count == 9:
      syntax = syntax + r'=*9:*'
      count = 1
    if i >= 1:
      syntax = (syntax + r'\t*')
    count = count + 1
    i = i + 1


  # Append ARB code to .ift filter for the purposes of extracting the sequences
  # from your custom databse.
  outFile.write('\n' + 'SEQUENCEAFTER   "*"'+'\n')
  outFile.write('SEQUENCESRT     ""'+'\n')
  outFile.write('SEQUENCEACI     "remove(".-0123456789 /")"'+'\n')
  outFile.write('SEQUENCECOLUMN 	0'+'\n')
  outFile.write('SEQUENCEEND     ">*"'+'\n\n')
  outFile.write('DONT_GEN_NAMES'+'\n\n')
  outFile.write('END             "//"')

  metaLabels.close()
  outFile.close()

if __name__ == "__main__":
  sys.exit(main())


