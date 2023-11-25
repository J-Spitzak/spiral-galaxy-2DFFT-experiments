#!/bin/bash

read -p "input file: " GALAXY 

cd ./$GALAXY

ls

read -p "input file: " FILE

ds9 $FILE &

read -p "x: " XPOS

read -p "y: " YPOS 

read -p "scale: " SCALE 

python3 ../image_cutout.py $FILE $XPOS $YPOS $SCALE $GALAXY.fits

python3 ../spiral_overlay.py $GALAXY.fits 

cp ../p2dfft ./

./p2dfft $GALAXY.fits

python3 ../p2pa $GALAXY.fits

more Results.csv

python3 ../spiral_overlay.py $GALAXY.fits

