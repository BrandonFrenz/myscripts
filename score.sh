#!/bin/sh
for i in *pdb
do
    ~/Desktop/Rosetta/main/source/bin/score_jd2.default.linuxgccrelease  \
        -s $i \
        -ignore_unrecognized_res \
        -edensity:mapfile $1 \
        -out:suffix $2 \
        -overwrite \
        -score:weights my.wts \
        -in::file::centroid_input\
        -mapreso 2 \


done
