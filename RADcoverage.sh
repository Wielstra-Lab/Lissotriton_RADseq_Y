#!/bin/bash

input=$1

loci_count=1

cat $input | while read locus base coverage
do
if [ $base = 10 ]; then
	
	while [ $loci_count -lt $locus ]; do
		echo -e  "$loci_count\t0"
		loci_count=$((loci_count+1))
	done

	echo -e "$locus\t$coverage"
	loci_count=$((loci_count+1))

fi
done
