#!/usr/bin/env python
###############################################################################
# Program: rename_tree_leaves.py
# Rev: 1.0
# Date: 8/20/12
# Author: Steve Essinger, Email: sde22@drexel.edu
#
# Description: Assigns unique ID from database to corresponding leaves on tree
#
# Example:
#
# rename_tree_leaves.py -r MFS_Aligin.fasta -i MFS_UID.fasta -l TreeLabels_orig.txt -o TreeLabels_Mapped_New.txt

import sys
import argparse

def main():

  if len(sys.argv) != 1:
    parser = argparse.ArgumentParser(description="Assigns unique ID from database to corresponding leaves on tree")
    parser.add_argument("-a", "--aligned-fasta", help="input aligned reference fasta file", required=True)
    parser.add_argument("-u", "--uid-fasta", help="input unique ID'd labeled fasta file", required=True)
    parser.add_argument("-l", "--labels", help="input tree labels file", required=True)
    parser.add_argument("-o", "--out", help="mapped output tree file", required=True)
    args = parser.parse_args()

    ref_fn  = args.aligned_fasta
    id_fn   = args.uid_fasta
    tree_fn = args.labels
    out_fn  = args.out
  else:
    ref_fn  = raw_input('Enter aligned fasta file:')
    id_fn   = raw_input('Enter UID fasta file:')
    tree_fn = raw_input('Enter input tree labels file:')
    out_fn  = raw_input('Enter output tree file:')

  refFile = open(ref_fn, 'r')
  idFile  = open(id_fn,  'r')
  giFile  = open(tree_fn,'r')
  outFile = open(out_fn, 'w')

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

if __name__ == "__main__":
  sys.exit(main())
