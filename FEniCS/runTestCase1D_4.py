""" Run this file."""
import numpy as np
import testCases as testCases
from testCaseWrapper import wrapper1D

if __name__ == '__main__':
    # select the parent folder as working directory !!!
    ModelOrder = [1, 3]
    nGridPoints = 641
    g = [0.0, 0.01, 0.3]
    for gi in g:
        tmp = 'testCase1D_4_g%1.3f'%gi
        tmp = tmp.replace('.', ',')
        solutionFolder = '../build/%s/FEniCS/'%tmp
        boundaryFilePrefix = '../build/%s/systemBoundary_%s'%(tmp, tmp)
        par = testCases.loadTestCase1D_4(gi)
        for N in ModelOrder:
            print('N: %d, g: %d'%(N, gi))
            u = wrapper1D(par, nGridPoints, N, boundaryFilePrefix, solutionFolder)
