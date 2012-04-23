#! /bin/bash

iters=100
inpath="../testdata/"
outpath="../analysis/"
files="114bp-5sp-test.phy 60bp-30sp-test.phy 300bp-5sp-test.phy 300bp-30sp-test.phy"

for f in $files; do

	outfile=`echo $f | cut -f1 -d .`"-all-1.txt"
	for (( i=0; i<iters; i++ )); do
		./dnaml1 $inpath$f 2>> $outpath$outfile > /dev/null

	done

done

for f in $files; do

	outfile=`echo $f | cut -f1 -d .`"-all-8.txt"
	for (( i=0; i<iters; i++ )); do
		./dnaml8 $inpath$f 2>> $outpath$outfile > /dev/null

	done

done

