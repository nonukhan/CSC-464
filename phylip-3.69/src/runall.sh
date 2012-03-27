#! /bin/bash

iters=1000
inpath="../testdata/"
outpath="../analysis/"
files="114bp-5sp.txt, 60bp-30sp.txt, 300bp-5sp.txt, 300bp-30sp.txt"

for f in $files; do

	outfile=`echo $f | cut -f1 -d .`"-all.txt"
	for (( i=0; i<iters; i++ )); do

		./dnaml $path$f >> $outpath$outfile

	done

done
