#! /bin/bash

#	Use to generate test sqquences 
#	./seq-gen <#base-pairs> <#sequences>
#


bp=$1
sp=$2

outfile=$bp"bp-"$sp"sp-test.phy"

echo $2 $1 > $outfile

for (( j=0; j<$sp; j++ )); do

	currline="Sequence$j	"
	for (( i=0; i<$bp; i++ )); do
		num=$RANDOM
		(( num %= 4 ))
		case $num in
			0)
			currline=$currline"A"
			;;
			1)
			currline=$currline"T"
			;;
			2)
			currline=$currline"C"
			;;
			3)
			currline=$currline"G"
			;;
		esac
	done

	echo $currline >> $outfile
done
