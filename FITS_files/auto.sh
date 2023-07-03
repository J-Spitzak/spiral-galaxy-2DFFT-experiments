#!/bin/bash

read -p "input file: " FILE

cd ./$FILE

ls

read -p "input file: " OTHER

cp $OTHER $FILE.fits

python3 ../spiral_overlay.py $FILE.fits 

cp ../p2dfft ./

./p2dfft $FILE.fits

python3 ../p2pa $FILE.fits

more Results.csv

python3 ../spiral_overlay.py $FILE.fits

