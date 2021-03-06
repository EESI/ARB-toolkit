
# inputs
# 1) Aligned fasta
# 2) Phylogenic Tree
# 3) Metadata fasta file
# 4) metadata labels file.
# 5) output directory

# outputs
# IFT import filter for arb
# Fasta with Unique Id's
# Tree with Unique Id's

set -e 
align_fa=$1
tree=$2
metadata=$3
metalabels=$4
output_dir=$5

if (( $# < 5 )); then  
	echo "$0: aligned.fasta tree metadata metadata-label output"
	exit
else


echo 
echo "VARIABLES:"
echo align_fa=$1
echo tree=$2
echo metadata=$3
echo metalabels=$4
echo output_dir=$5
fi

# Make our output directory
if [[ ! -d "$output_dir" ]]; then
	mkdir -p "$output_dir"
fi
accession="$output_dir/accession"
UID_fa="$output_dir/$(basename $align_fa)"
UID_tree="$output_dir/$(basename $tree)"
treelabels_original="$output_dir/extracted_labels"
treelabels_mapped="$output_dir/mapped_labels"
ift_filter="$output_dir/$(basename $align_fa).ift"
echo


echo "getting accession numbers"
getAccession.py -i $metadata -o $accession # works
echo "adding Unique ID's to aligned fasta"
addUIDtoFasta.py -i $accession -a $align_fa -o $UID_fa # works
echo "extracting leaf names"
extract_leaf_names.py -i $tree -o $treelabels_original # doesn't work!
echo "building new tree leaf names"
rename_tree_leaves.py -a $align_fa -u $UID_fa -l $treelabels_original -o $treelabels_mapped
echo "converting tree to new names"
convert_leaf_names.py -i $tree -l $treelabels_mapped -o $UID_tree
echo "building ARB ift filter"
build_ift_from_metalabels.py -i $metalabels -o $ift_filter
