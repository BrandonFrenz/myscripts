#!/bin/bash
./scons.py $1 $2 $3 $4
play --no-show-progress --null --channels 1 synth 1 sine 400
