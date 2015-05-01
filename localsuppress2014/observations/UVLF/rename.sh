#!/bin/bash

files=Duncan*

for f in $files
do
	fname="${f%.*}"
	new=${fname}.txt
	mv $f $new
done
