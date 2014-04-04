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
    parser.add_argument("-o", "--out", help="mapped mapped tree labels file", required=True)
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

  labels_fh= open(tree_fn,'r')
  out_fh   = open(out_fn, 'w')

  headers = []
  ids = []

  # add all our headers to headers
  with open(ref_fn, 'r') as fh:
    for line in fh:
      if line.startswith('>'):
        headers.append(line.strip())

  # Add our Uids's to ids
  with open(id_fn, 'r') as fh:
    for line in fh:
      if line.startswith('>'):
        # remove > and \n
        ids.append(line[1:-1])

  # iterate over each line in our tree labels file
  for line in labels_fh:
    orig = line.strip()
    orig = orig.partition('|')
    line = line.partition('_')
    temp = line[2].partition('|')
    line = temp[0]

    for i in range(0,len(headers)):
      if line.strip() in headers[i]:
        if temp[1]:
          out_fh.write(str(orig[0])+'A'+'\t'+str(ids[i].strip())+'\n')
        else:
          out_fh.write(str(orig[0])+'\t'+str(ids[i].strip())+'\n')
        continue

    #if i == len(REF)-1:
    #  print orig[0].strip()

  out_fh.close()

if __name__ == "__main__":
  sys.exit(main())
