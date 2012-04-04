#! /bin/bash

if [ $1 == "cu" ]; then
	files="dnaml.c phylip.c seq.c ../../simple-timing/timekeeper.c"
else 
	if [ $1 == "c" ]; then
		files="dnaml.cu phylip.cu seq.cu ../../simple-timing/timekeeper.cu"

	fi
fi

for f in $files; do
	
	name=`echo $f | cut -f1 -d c`$1
	mv $f $name

done
