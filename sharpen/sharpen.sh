#!/bin/bash
~/Desktop/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
    -parser:protocol bfactor.xml \
    -parser:script_vars map=$1 \
    -nstruct 1 \
    -in::file::fasta $2
