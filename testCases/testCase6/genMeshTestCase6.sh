#!/bin/bash
# run this script with your command line knowing gmsh and dolfin-convert (e.g. through conda activate ...)

# refine by splitting
 cd ../../tools
./genMeshRefine.sh -r 1 -s ../testCases/testCase6/ -n meshTestCase6 -d ../build/testCase6/;
 cd ../testCases/testCase6/

mv meshTestCase6Normals.mat ../../build/testCase6/


