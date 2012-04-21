#! /bin/bash

iters=1000
inpath="../testdata/"
outpath="../analysis/"
files="114bp-5sp-test.phy 60bp-30sp-test.phy 300bp-5sp-test.phy 300bp-30sp-test.phy"

for f in $files; do

	outfile=`echo $f | cut -f1 -d .`"-all.txt"
	for (( i=0; i<iters; i++ )); do
		./dnaml $inpath$f 2>> $outpath$outfile > /dev/null

	done

done
