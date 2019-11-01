""" Run this file to compute the solution of the second-order formulation
of the PN equations for test case 1 (approximated by the finite element method).

%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
"""


import numpy as np
import testCases as testCases
from testCaseWrapper import wrapper1D

if __name__ == '__main__':
    solutionFolder = '../build/testCase1/FEniCS/'
    boundaryFilePrefix = '../build/testCase1/systemBoundary_testCase1'
    ModelOrder = np.arange(1, 22, 2, dtype='int')
    nGridPoints = 1 + np.array([10, 20, 40, 80, 160, 320, 640])
    par = testCases.loadTestCase1()
    for n in nGridPoints:
        for N in ModelOrder: 
            print('N: %d, grid: %d'%(N, n))
            u = wrapper1D(par, n, N, boundaryFilePrefix, solutionFolder)
