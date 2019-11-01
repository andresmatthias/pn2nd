"""Run this file."""
import testCases as testCases
from testCaseWrapper import wrapper2D
    

if __name__ == '__main__':
    # select the parent folder as working directory !!!
    solutionFolder = '../build/testCase2D_4/FEniCS/'
    boundaryFilePrefix = '../build/testCase2D_4/systemBoundary_testCase2D_4'
    ModelOrder = [1, 3, 5, 7]
    meshNames = ['../build/testCase2D_4/meshTestCase2D_4_0',
                 '../build/testCase2D_4/meshTestCase2D_4_1']
    par = testCases.loadTestCase2D_4()
    for m in meshNames:
        for N in ModelOrder: 
            u = wrapper2D(par, m, N, boundaryFilePrefix, solutionFolder)
