#!/bin/bash

# execute this file within this folder

refinesteps=1
destination=../build/testCase2D_1/
sourcefolder=../testCases/testCase2D_1/
meshname_prefix=meshTestCase2D_1

# parse arguments
while getopts r:s:n:d: option
do
case "${option}"
in
r) refinesteps=${OPTARG};;
s) sourcefolder=${OPTARG};;
n) meshname_prefix=${OPTARG};;
d) destination=${OPTARG};;
esac
done


# create initial triangulation
meshname0_prefix=$meshname_prefix$(echo "_0")
lastmeshname=$sourcefolder$meshname0_prefix$(echo ".msh")
echo $sourcefolder$meshname_prefix.geo 
gmsh $sourcefolder$meshname_prefix.geo -2 -o $lastmeshname
#gmsh $sourcefolder$meshname_prefix.geo -3 -o $lastmeshname

meshnamelist=($meshname0_prefix)

# refine
for (( c=1; c<=$refinesteps; c++ ))
do
    refinedmeshname_prefix=$meshname_prefix$(echo "_")$c
    refinedmeshname=$sourcefolder$refinedmeshname_prefix$(echo ".msh")
    gmsh $lastmeshname -refine -o $refinedmeshname
    lastmeshname=$refinedmeshname
    meshnamelist+=($refinedmeshname_prefix)
done

meshnamelist_length=${#meshnamelist[@]}

# move and convert
for (( i=0; i<${meshnamelist_length}; i++ ));
do
    meshname_prefix=${meshnamelist[$i]}
    meshname=$sourcefolder$meshname_prefix$(echo ".msh")
    movedmeshname_prefix=$destination$meshname_prefix
    movedmeshname=$movedmeshname_prefix$(echo ".msh")
    mv $meshname $movedmeshname
    dolfin-convert $movedmeshname $movedmeshname_prefix$(echo ".xml")
done

