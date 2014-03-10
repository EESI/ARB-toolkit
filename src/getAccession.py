#!/usr/bin/env python
###############################################################################
# Program: getAccession.py
# Rev: 1.0
# Date: 8/2/12
# Author: Steve Essinger, Email: sde22@drexel.edu
#
# Description: Maps unique ID from database to each alignment
###############################################################################

import sys
import argparse

def main():

  if len(sys.argv) != 1:
    parser = argparse.ArgumentParser(description="Maps unique metadata id's to accession list")
    parser.add_argument("-i", "--input-metadata", help="input metadata file", required=True)
    parser.add_argument("-o", "--output-accession", help="output accession list file", required=True)
    args = parser.parse_args()

    inFile = args.input_metadata
    outFile = args.output_accession
  else:
    inFile = raw_input('Enter input metadata file:')
    outFile   = raw_input('Enter output accession list file:')

  inFile = open(inFile,'r')
  outFile = open(outFile,'w')

  for line in inFile:
    if line.startswith('>'):
      temp = line.split('\t',12)
      revised = temp[11].partition('_')
      if revised[1]:
        revised = revised[0]+'|A'
      else:
        revised = revised[0]
      outFile.write(str(temp[1])+'\t'+str(revised)+'\n')
  inFile.close()
  outFile.close()

if __name__ == "__main__":
  sys.exit(main())
