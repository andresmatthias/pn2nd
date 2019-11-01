""" Run this file to compute the solution of the second-order formulation
of the PN equations for test case 5 (approximated by the finite element method).

%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
"""

import testCases as testCases
from testCaseWrapper import wrapper2D
    

if __name__ == '__main__':
    ModelOrder = [1, 3, 5, 7]
    meshNames = ['../build/testCase5/meshTestCase5_0',
                 '../build/testCase5/meshTestCase5_1']
    g = [0.0, 0.5]
    for gi in g:
        tmp = 'testCase5_g%1.3f'%gi
        tmp = tmp.replace('.', ',')
        solutionFolder = '../build/%s/FEniCS/'%tmp
        boundaryFilePrefix = '../build/testCase5/systemBoundary_testCase5'
        par = testCases.loadTestCase5(gi)
        for m in meshNames:
            for N in ModelOrder: 
                u = wrapper2D(par, m, N, boundaryFilePrefix, solutionFolder)
