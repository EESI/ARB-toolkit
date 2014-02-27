#!/usr/bin/env python 

# *** Gail edit -- makes sure names_scientific.dmp only contain scientific names***
# print """ 
# Need file with taxids as input argument
#This program reads the NCBI taxonomy database (text file format) and
#reconstructs the whole taxonomy tree of life. Once loaded, the program
#allows you to look up for the specific subtree of a given taxid, or
#save the whole tree using the newick format.

#Note that reconstructing the from the raw NCBI format may take some
#minutes.
#"""

#11/30/11 S.Essinger
#This is a modified version of Gail's code. It is designed to be used as a
#python module. There are two functions:
#build_tax_tree() reconstructs the whole taxonomy tree of life. it required as
#preprocessing step for the get_taxonomy() function
#get_taxonomy() input: taxid, id2node, id2name; output: returns the rank and
#names of the taxid's taxonomy.

'''
9.19.2012 E.Reichenberger
Re-downloading the names.dmp from NCBI will have scientific and non-scientific names within the file. 
To ensure you have only the scientific names:
cat names.dmp | grep "scientific name" > names_scientific.dmp

'''


def build_tax_tree():
        import os
        import sys 
        from string import strip
        from ete2 import TreeNode, Tree
        #print sys.argv[1]
        #if len(sys.argv) == 1:
        #	print "Usage:  taxid2lineage file_with_taxids.txt"
        #else:
        #	f = open(sys.argv[1], 'r')

        # This sets Unbuffered stdout/auto-flush
        sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

        id2node= {}
        id2rank={}
        node2parentid = {}
        all_ids = set([])
        all_nodes = []
        id2name= {}

        # Loads info from NCBI taxonomy files
        if os.path.exists("nodes.dmp"):
            NODESFILE = open('nodes.dmp')
        elif os.path.exists("nodes.dmp.bz2"):
            import bz2
            NODESFILE = bz2.BZ2File('nodes.dmp.bz2')
        else:
            print '"nodes.dmp" file is missing. Try to downloaded from: '

        if os.path.exists("names_scientific.dmp"):
            NAMESFILE = open('names_scientific.dmp')
        elif os.path.exists("names_scientific.dmp.bz2"):
            import bz2
            NAMESFILE = bz2.BZ2File('names_scientific.dmp.bz2')
        else:
            print '"names_scientific.dmp" file is missing. Try to downloaded from: '

        # Reads taxid/names transaltion
        #print 'Loading species names from "names_scientific.dmp" file...',
        for line in NAMESFILE:
            line = line.strip()
            fields = map(strip, line.split("|"))
            nodeid, name = fields[0], fields[1]
            id2name[nodeid] = name


        # Reads node connections in nodes.dmp
        #print 'Loading node connections form "nodes.dmp" file...', 
        for line in NODESFILE:
            line = line.strip()
            fields = map(strip, line.split("|"))
            nodeid, parentid,rankid = fields[0], fields[1], fields[2]
            id2rank[nodeid]=rankid
            if nodeid =="" or parentid == "":
                raw_input("Wrong nodeid!")

            # Stores node connections
            all_ids.update([nodeid, parentid])

            # Creates a new TreeNode instance for each new node in file
            n = TreeNode()
            # Sets some TreeNode attributes
            n.add_feature("name", id2name[nodeid])
            n.add_feature("taxid", nodeid)
            n.add_feature("rank",id2rank[nodeid])

            # updates node list and connections
            node2parentid[n]=parentid
            id2node[nodeid] = n
        #print len(id2node)

        # Reconstruct tree topology from previously stored tree connections
        #print 'Reconstructing tree topology...'
        for node in id2node.itervalues():
            parentid = node2parentid[node]
            parent = id2node[parentid]
            if node.taxid == "1":
                t = node
            else:
                parent.add_child(node)

        return id2node, id2name

def get_taxonomy(taxon_id,id2node,id2name):

        # Let's play with the tree
        def get_track(node):
            ''' Returns the taxonomy track from leaf to root'''
            track = []
            itsrank=[]
            while node is not None:
                if node.rank != "no rank":
                        track.append(node.name) # You can add name or taxid
                        itsrank.append(node.rank)
                node = node.up
            return track,itsrank  

        taxid = None
        track = None
        itsrank = None
        taxid = str(taxon_id)
        if taxid in id2name:
                target_node = id2node[taxid]
                track,itsrank=get_track(target_node)
	#track.reverse()
	#itsrank.reverse()
        return track, itsrank
