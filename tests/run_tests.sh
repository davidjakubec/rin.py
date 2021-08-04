#!/bin/bash

for structure_file in *.cif
do
	echo "Processing file:" ${structure_file}
	../rin.py ${structure_file}
	../rin.py ${structure_file} --node_atom_selection CB
	../rin.py ${structure_file} --node_atom_selection nonhydrogen --cutoff 4.5
	../rin.py ${structure_file} --node_atom_selection nonhydrogen --include_hetatms --cutoff 4.5
done
