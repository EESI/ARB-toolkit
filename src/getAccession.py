#!/usr/bin/env python
###############################################################################
# Program: getAccession.py
# Rev: 1.0
# Date: 8/2/12
# Author: Steve Essinger, Email: sde22@drexel.edu
#
# Description: Maps unique ID from database to each alignment
###############################################################################

inFile = open('MFS_metaData.txt','r')
outFile = open('ListAccessions.txt','w')

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

inACC = open('ListAccessions.txt','r')
inALN = open('MFS_Align.fasta','r')
outALN = open('MFS_UID.fasta','w')

UID = []
ACC = []
for line in inACC:
  line = line.strip()
  line = line.partition('\t')
  hack = line[2].partition('|')
  UID.append(str(line[0]))
  ACC.append(str(hack[0]))
inACC.close()

for line in inALN:
  if line.startswith('>'):
    outALN.write('\n')
    temp = line.partition('pdb|')
    if len(temp[1]) == 0:
      temp = line.partition('sp|')
    temp = temp[2].partition('|')
    temp = temp[0].partition('.')
    try:
      key = ACC.index(temp[0])
      plug = UID[key]
      outALN.write('>'+str(plug)+'\n')
    except ValueError:
      key = 'MISSING_IN_DATABASE'
      outALN.write('>'+str(key)+'\n')
  else:
    outALN.write(line)
inALN.close()
outALN.close()
