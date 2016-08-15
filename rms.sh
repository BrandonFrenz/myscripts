#!/bin/sh
for i in *pdb
do
    ~/Desktop/Rosetta/main/source/bin/score_jd2.default.linuxgccrelease  \
        -in:file:native $1 \
        -s $i \
        -ignore_unrecognized_res \
        -edensity:mapfile $2 \
        -fastdens_wt 20 \
        -out:suffix $3 
done
