#!/bin/sh
for i in *pdb
do
    perl ~/Desktop/scripts/batch_select_chains.pl AB $i
done
