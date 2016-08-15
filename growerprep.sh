#!/bin/bash
/home/brandon/Desktop/Rosetta/main/source/bin/grower_prep.default.linuxgccrelease \
    -in::file::fasta $2 \
    -pdb $1 \
    -fragsizes 3 9 \
    -fragamounts 100 20 \
