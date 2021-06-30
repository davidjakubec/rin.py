#!/bin/bash

for structure_file in *.cif
do
	echo "Processing file:" ${structure_file}
	../rin.py ${structure_file}
done
