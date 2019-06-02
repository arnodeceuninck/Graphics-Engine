#!/bin/sh

# Author: Arno Deceuninck
# Year: 2019

mkdir test_output 2>/dev/null
rm -r test_output/*  # cleanup the old files

##################################################

cd l_systems

for i in *.ini; do
    echo Testing $i
    ../engine $i || break
done

mv *.bmp ../test_output
cp *.png ../test_output
# compare image1 image2 -compose src diff.png
cd ../

###################################################

cd line_drawings

for i in *.ini; do
    echo Testing $i
    ../engine $i || break
done

mv *.bmp ../test_output
cp *.png ../test_output
cd ../

###################################################

cd wireframes

for i in *.ini; do
    echo Testing $i
    ../engine $i
done

mv *.bmp ../test_output
cp *.png ../test_output
cd ../

###################################################

cd z_buffered_wireframes

for i in *.ini; do
    echo Testing $i
    ../engine $i
done

mv *.bmp ../test_output
cp *.png ../test_output
cd ../

###################################################

cd z_buffering

for i in *.ini; do
    echo Testing $i
    ../engine $i
done

mv *.bmp ../test_output
cp *.png ../test_output
cd ../

###################################################

cd 3d_fractals

for i in *.ini; do
    echo Testing $i 
    ../engine $i 
done

mv *.bmp ../test_output
cp *.png ../test_output 
cd ../

###################################################

cd ambient_light

for i in *.ini; do
    echo Testing $i 
    ../engine $i 
done

mv *.bmp ../test_output
cp *.png ../test_output 
cd ../

###################################################

cd shadowing

for i in *.ini; do
    echo Testing $i 
    ../engine $i 
done

mv *.bmp ../test_output
cp *.png ../test_output 
cd ../

###################################################

cd spheres_and_cylinders

for i in *.ini; do
    echo Testing $i 
    ../engine $i 
done

mv *.bmp ../test_output
cp *.png ../test_output 
cd ../

###################################################



cd test_output

bmp=".bmp"
png=".png"
compared=".diff.png"

for i in *.bmp; do
    filename="${i##*/}"
    filename="${filename%.*}"
    filenamepng="$filename$png"
    filenamebmp="$filename$bmp"
    filenamecompared="$filename$compared"
    echo $filenamepng comparing with $filenamebmp
    compare $filenamepng $filenamebmp -compose src $filenamecompared
done

mkdir compared
mv *.diff.png compared

mkdir empty

cd ../

# References
# https://stackoverflow.com/questions/5132749/diff-an-image-using-imagemagick
# https://stackoverflow.com/questions/14505047/loop-through-all-the-files-with-a-specific-extension
