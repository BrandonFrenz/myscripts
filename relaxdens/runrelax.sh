#!/usr/bin/bash
/Users/brandonfrenz/rosetta-main/source/bin/rosetta_scripts.default.macosclangrelease\
    -database /Users/brandonfrenz/rosetta-main/database/\
    -s $1\
    -edensity:mapfile $2\
    -parser:protocol relax_into_density.xml\
    @relax_into_density.flags
