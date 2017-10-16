#!/bin/sh
for i in *pdb
do
    ~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease  \
    -s $i \
    -parser:protocol $1 \
    -edensity:mapfile $2\
    -default_max_cycles 200 \
    -ignore_unrecognized_res \
    -nstruct 1 \
    -overwrite \
    -missing_density_to_jump \
    -beta

done
