#!/bin/sh

~/Desktop/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
    -parser:protocol $1 \
    -edensity:mapfile $3\
    -s $2 \
    -default_max_cycles 200 \
    -ignore_unrecognized_res \
    -nstruct 1 \
    -overwrite \
    -parser::script_vars readbeams=$5 \
    -parser::script_vars beams=$6 \
    -parser::script_vars steps=$7 \
    -parser::script_vars pcount=$8 \
    -parser::script_vars filterprevious=$9 \
    -parser::script_vars filterbeams=${10} \
    -in::file::native $4 \
    -mapreso 2 \
    -beta \
    -missing_density_to_jump\
    #-cenrot \
    #-edensity:centroid_density_mass 3 \

