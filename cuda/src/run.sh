#! /bin/bash

iters=1000
path="../testdata/"


for (( i=0; i<iters; i++ )); do

	./dnaml $1 

done


