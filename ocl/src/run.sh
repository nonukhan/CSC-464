#! /bin/bash

iters=1000
path="../testdata/"
files="114bp-5sp.txt, 60bp-30sp.txt, 300bp-5sp.txt, 300bp-30sp.txt"

#for f in $files; do

	for (( i=0; i<iters; i++ )); do

		./dnaml $1 

	done

#done
