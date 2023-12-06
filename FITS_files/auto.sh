#!/bin/bash

read -p "input galaxy: " GALAXY 

cd ./$GALAXY

ls
echo ""
echo "your input file should look something like this:"
echo "fc_184.706682+14.416509_sdss(dr7)_g.fits"
echo "remember the last letter"
echo ""

read -p "input file: " FILE
read -p "last letter: " BAND
echo ""

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

