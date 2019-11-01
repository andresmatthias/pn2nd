#!/bin/bash
# run this script with your command line knowing gmsh and dolfin-convert (e.g. through conda activate ...)

# refine by splitting
 cd ../../tools
./genMeshRefine.sh -r 1 -s ../testCases/testCase4/ -n meshTestCase4 -d ../build/testCase4/;
 cd ../testCases/testCase4/

# refine by halving characteristic length
#python3 ../../tools/genMesh.py ./ meshTestCase4 -d ../../build/testCase4/ -l 3;


