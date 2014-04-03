#!/usr/bin/env python
#
#  Program: extract_leaf_names.py
#
#  Description: extracts clade names from a newick phylogenic 
# 
#  Example: 
#
#  extract_leaf_names.py -i mfs.tree -o treelabels.txt 
import Bio.Phylo as p
import sys
import argparse

def main():

	parser = argparse.ArgumentParser(description="extracts clade names from a newick phylogenic")
	parser.add_argument("-i", "--input", help="input tree file", required=True)
	parser.add_argument("-o", "--output", help="extracted clade names output file", required=True)
	args = parser.parse_args()

	tree = p.read(args.input, 'newick')
	out_fh = open(args.output, 'w')

	for clade in tree.find_elements({}):
		if clade.name is not None:
			out_fh.write(str(clade.name) + "\n")

if __name__ == "__main__":
	sys.exit(main())
