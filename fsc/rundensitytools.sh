#!/usr/bin/bash
/Users/brandonfrenz/rosetta-main/source/bin/density_tools.default.macosclangrelease\
    -s $1 \
    -edensity::mapfile $3 \
    -edensity::alt_mapfile $4 \
    -edensity::cryoem_scatterers true \
    -denstools:lowres 20 -edensity:mapreso $2 -nresbins 100 -denstools:verbose -denstools:perres
