#!/bin/bash

echo "give tsv list?"
read filename

while read -r line
do
    tfile=$line
    echo "tsv file is: $tfile"

	cat $tfile | sed 's/\,/\./g' > somecopy.txt 
	/bin/mv somecopy.txt $tfile

done < "$filename"







