#!/bin/bash
# run this script with your command line knowing gmsh and dolfin-convert (e.g. through conda activate ...)

# refine by splitting
 cd ../../tools
./genMeshRefine.sh -r 1 -s ../testCases/testCase5/ -n meshTestCase5 -d ../build/testCase5/;
 cd ../testCases/testCase5/

# refine by halving characteristic length
#python3 ../../tools/genMesh.py ./ meshTestCase5 -d ../../build/testCase5/ -l 3;


