#!/bin/bash
check () { 
	if [[ "$1" -ne "0" ]]; then
		echo "FAILED $3"
    if [[ $debug ]]; then 
			less <<< $2
			exit 1
		else 
			exit 1 
		fi
	else
		echo "SUCCESS $3"
	fi
}


out="test"

echo "RUNNING PIPELINE"
../src/arb_import_pipeline.sh ../example/MFS_Align.fasta ../example/MFS.tree ../example/MFS_metaData.txt ../example/MFS_Field_labels.txt $out 
if [[ "$?" -ne "0" ]]; then
	echo "PIPELINE FAILED"
fi


# tests
res=$(diff ../example/ListAccessions.txt $out/accession)
check $? "$res" "checking accession number "

res=$(diff ../example/MFS_UID.fasta $out/MFS_Align.fasta)
check $? "$res" "checking UID fasta file "

res=$( diff -w <( sort ../example/TreeLabels_Orig.txt ) <( sort $out/extracted_labels ))
check $? "$res" "checking extracted tree labels "

res=$( diff -w <( sort ../example/TreeLabels_Mapped_New.txt ) <( sort $out/mapped_labels ))
check $? "$res" "checking mapped labels "

res=$( diff -w ../example/MFS.tree $out/MFS.tree )
check $? "$res" "checking tree file"

res=$( diff -w ../example/MFS_import_filter.ift $out/MFS_Align.fasta.ift )
check $? "$res" "checking ARB ift filter"
