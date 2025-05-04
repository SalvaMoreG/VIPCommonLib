#!/bin/bash

echo "give tsv list?"
read filename

cat $filename | sed 's/\,/\./g' > somecopy.txt 
/bin/mv somecopy.txt $filename








