#!/bin/bash
check () { 
	if [[ "$1" -ne "0" ]]; then
		echo -n "FAILED "
    # less <<< $2
	else
		echo -n "SUCCESS "
	fi
}


out="test"
mkdir  $out

echo "RUNNING PIPELINE"
../src/arb_import_pipeline.sh ../example/MFS_Align.fasta ../example/MFS.tree ../example/MFS_metaData.txt ../example/MFS_Field_labels.txt $out 


# tests
res=$(diff ../example/ListAccessions.txt $out/accession)
check $? $res
echo "checking accession number "

res=$(diff ../example/MFS_UID.fasta $out/MFS_Align.fasta)
check $? $res
echo -n "checking UID fasta file "

res=$( diff -w <( sort ../example/TreeLabels_Orig.txt ) <( sort $out/extracted_labels ))
check $? $res
echo "checking extracted tree labels "

res=$( diff -w <( sort ../example/TreeLabels_Mapped_New.txt ) <( sort $out/mapped_labels ))
check $? $res
echo "checking mapped labels "

res=$( diff -w ../example/MFS.tree $out/MFS.tree )
check $? $res
echo "checking tree file"

res=$( diff -w ../example/MFS_import_filter.ift $out/MFS_Align.fasta.ift )
check $? $res
echo "checking ARB ift filter"
