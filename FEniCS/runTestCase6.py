""" Run this file to compute the solution of the second-order formulation
of the PN equations for test case 6 (approximated by the finite element method).

%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
"""

import testCases as testCases
from testCaseWrapper import wrapper2D
    
if __name__ == '__main__':
    solutionFolder = '../build/testCase6/FEniCS/'
    boundaryFilePrefix = '../build/testCase6/systemBoundary_testCase6'
    ModelOrder = [1, 3, 5, 7]
    meshNames = ['../build/testCase6/meshTestCase6_0',
                 '../build/testCase6/meshTestCase6_1']
    par = testCases.loadTestCase6() 
    for m in meshNames:
        for N in ModelOrder: 
            u = wrapper2D(par, m, N, boundaryFilePrefix, solutionFolder)
