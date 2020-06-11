#!/usr/bin/bash
/home/jovyan/main/source/bin/rosetta_scripts.default.linuxgccrelease\
    -database /home/jovyan/main/database/\
    @rayfineflags\
    -parser:protocol targeted.xml\
    -edensity:mapfile croppedfull.mrc\
    -out:suffix $2\
    -s $1
