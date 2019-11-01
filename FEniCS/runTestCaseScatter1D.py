""" Run this file to compute the solution of the second-order formulation
of the PN equations for test case 1 (approximated by the FEM method).

%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
"""


import numpy as np
import testCases as testCases
from testCaseWrapper import wrapper1D

if __name__ == '__main__':
    solutionFolder = '../build/testCaseScatter/FEniCS/'
    boundaryFilePrefix = '../build/testCaseScatter/systemBoundary_testCaseScatter'
    ModelOrder = np.array([1, 3, 5])
    nGridPoints = 1 + np.array([100])
    g = 0.4
    par = testCases.loadTestCaseScatter1D(g)
    for n in nGridPoints:
        for N in ModelOrder: 
            print('N: %d, grid: %d'%(N, n))
            u = wrapper1D(par, n, N, boundaryFilePrefix, solutionFolder)
