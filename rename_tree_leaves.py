#!/usr/bin/env python
###############################################################################
# Program: rename_tree_leaves.py                                                    
# Rev: 1.0                                                                    
# Date: 8/20/12                                                                
# Author: Steve Essinger, Email: sde22@drexel.edu
#
# Description: Assigns unique ID from database to corresponding leaves on tree
###############################################################################

refFile = open('MFS_Align.fasta','r')
idFile = open('MFS_UID.fasta','r')
giFile = open('TreeLabels_Orig.txt','r')
outFile = open('TreeLabels_Mapped_New.txt','w')

REF = []
for line in refFile:
    if line.startswith('>'):
        REF.append(line.strip())
refFile.close()

ID = []
for line in idFile:
    if line.startswith('>'):
        temp = line.partition('>')
        ID.append(temp[2].strip())
idFile.close()

for line in giFile:
    orig = line.strip()
    orig = orig.partition('|')
    line = line.partition('_')
    temp = line[2].partition('|')
    line = temp[0]
    for i in range(0,len(REF)):
        if line.strip() in REF[i]:
            if temp[1]:
                outFile.write(str(orig[0])+'A'+'\t'+str(ID[i].strip())+'\n')
            else:
                outFile.write(str(orig[0])+'\t'+str(ID[i].strip())+'\n')
            continue
    if i == len(REF)-1:
        print orig[0].strip()

giFile.close()
outFile.close()
