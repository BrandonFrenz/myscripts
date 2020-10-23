#!/usr/bin/bash
$ROSETTA/source/bin/rosetta_scripts.default.macosclangrelease\
    -database $ROSETTA/database/\
    -s $1\
    -edensity:mapfile $2\
    -parser:protocol relax_into_density.xml\
    @relax_into_density.flags
