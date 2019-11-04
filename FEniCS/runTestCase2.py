""" Run this file to compute the solution of the second-order formulation
of the PN equations for test case 2 (approximated by the finite element method).

%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%
"""


import numpy as np
import testCases as testCases
from testCaseWrapper import wrapper1D

if __name__ == '__main__':
    solutionFolder = '../build/testCase2/FEniCS/'
    boundaryFilePrefix = '../build/testCase2/systemBoundary_testCase2'
    ModelOrder = np.arange(1, 22, 2, dtype='int')
    nGridPoints = 1 + np.array([10, 20, 40, 80, 160, 320, 640])
    par = testCases.loadTestCase2()
    for n in nGridPoints:
        for N in ModelOrder: 
            print('N: %d, grid: %d'%(N, n))
            u = wrapper1D(par, n, N, boundaryFilePrefix, solutionFolder)
