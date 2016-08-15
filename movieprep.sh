#!/bin/bash

for filename in *pdb; do
    z=$((z+1));
    mv $filename testing$z".pdb";
done

