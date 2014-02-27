A tutorial using ARB\_Toolkit is available:
http://www.ece.drexel.edu/gailr/EESI/tutorial.php

### Description: ###
Researchers are perpetually amassing biological sequence data. The
computational approaches employed by ecologists for organizing this data (e.g.
alignment, phylogeny, etc.) typically scale nonlinearly in execution time with
the size of the dataset. This often serves as a bottleneck for processing
experimental data since many molecular studies are characterized by massive
datasets. To keep up with experimental data demands, ecologists are forced to
choose between continually upgrading expensive in-house computer hardware or
outsourcing the most demanding computations to the cloud. Outsourcing is
attractive since it is the least expensive option, but does not necessarily
allow direct user interaction with the data for exploratory analysis. Desktop
analytical tools such as ARB are indispensable for this purpose, but they do
not necessarily offer a convenient solution for the coordination and
integration of datasets between local and outsourced destinations. Therefore,
researchers are currently left with an undesirable tradeoff between
computational throughput and analytical capability. To mitigate this tradeoff
we introduce a software package to leverage the utility of the interactive
exploratory tools offered by ARB with the computational throughput of
cloud-based resources. Our pipeline serves as middleware between the desktop
and the cloud allowing researchers to form local custom databases containing
sequences and metadata from multiple resources and a method for linking data
outsourced for computation back to the local database. A tutorial
implementation of the toolkit is provided in the supplementary material.

### Tools: ###

example usage:

		build_ift_from_metalabels.py -i ../arb/MFS_Field_labels.txt -o mfs-importer.ift

		rename_tree_leaves.py -r MFS_Align.fasta -i MFS_UID.fasta -l TreeLables_orig.txt -o treelabels_mapped.txt


## Usage: ##

### generate ift ###

first generate your ift from a file called Field labels. this is a file that
contains the label of each field in your database each on it's own line. If
you look in examples there is an example of this.

		build_ift_from_metalabels.py -i ../arb/MFS_Field_labels.txt -o mfs-importer.ift
		
now you need to copy your import filter to the $ARBHOME/lib/import folder



### Contributors ###
Steve Essinger
Erin Reichenberger
Chris  Blackwood
Gail Rosen
Calvin Morrison



