"""Solve the second-order formulation of the PN equations.

For details on how to use the defined functions, see
    ./testCaseWrapper.py
    ../unitTests.py

%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%
"""


import fenics as fe
import importlib
import numpy as np
from boundarySourceOddMoments import computeBoundarySourceOddMoments
import scipy.io as sio
import sys
sys.path.append('../build/PN2ndModels/')


def outerNormalVectors3DMat2List(outerNormalsAsMat):
    """Convert tensor of outer normal vectors (export from Matlab) to python list."""
    outerNormalsAsList = []
    for k in range(outerNormalsAsMat.shape[1]):
        outerNormalsAsList.append(outerNormalsAsMat[:, k].astype('float'))
    
    return outerNormalsAsList
    

def getNumberOfEvenBasisFunctions(N, spatialDimension):
    """Get number of even basis functions (real spherical harmonic S_l^m of
    degree l and order m even <=> l even) for given model order N.
    
    Depending on spatialDimension neglect certain basis functions due to
    symmetry assumptions.
    """
    if spatialDimension == 1:
        nEvenMoments = (N + 1) / 2
    elif spatialDimension == 2:
        nEvenMoments = 0
        for k in range(1, N + 1, 2):
            nEvenMoments += k
    elif spatialDimension == 3:
        nEvenMoments = 0
        for k in range(1, N + 1, 2):
            nEvenMoments += 2 * k - 1
    else:
        raise ValueError('Invalid spatial dimension!')

    return int(nEvenMoments)


def solvePN2nd(par, N_PN, mesh, boundaryFilePrefix, ds=fe.ds):
    """Assemble and solve the second-order formulation of the PN equations."""
    def initWeakFormulation():
        """Define function spaces etc. to initialize weak formulation."""
        P = fe.FiniteElement('P', mesh.ufl_cell(), 1)
        V = fe.FunctionSpace(mesh, 'P', 1)
        auxP = []
        for k in range(nEvenMoments):
            auxP.append(P)
        VVec = fe.FunctionSpace(mesh, fe.MixedElement(auxP))
        auxV = []
        for k in range(nEvenMoments):
            auxV.append(V)
        assigner = fe.FunctionAssigner(auxV, VVec)

        v = fe.TestFunctions(VVec)
        u = fe.TrialFunctions(VVec)
        uSol = fe.Function(VVec)
        
        uSolComponents = []
        for k in range(nEvenMoments):
            uSolComponents.append(fe.Function(V))
        
        # FEniCS work-around 
        if N_PN == 1: 
            solAssignMod = lambda uSolComponents, uSol: [assigner.assign(uSolComponents[0], uSol)]
        else:
           solAssignMod = lambda uSolComponents, uSol: assigner.assign(uSolComponents, uSol)
        
        return u, v, V, uSol, uSolComponents, solAssignMod


    def defineWeakSystemDomain():
        """Define system of equations of the weak formulation of 
        the second-orderformulation of the PN equations on the domain .
        """
        sMD = importlib.import_module(
            "systemMatricesDomain_d_%d_k_%s_N_%d" % (par['spatialDimension'], par['kernelName'], N_PN))
        systemMatricesDomain = sMD.getSystemMatricesDomain(par)

        lhs = 0.0 * v[0] * fe.dx(V)
        for i in range(nEvenMoments):
            for j in range(nEvenMoments):
                if (par['spatialDimension'] == 1):
                    lhs += systemMatricesDomain['Kzz'][i][j] * fe.dot(fe.grad(u[j])[0], fe.grad(v[i])[0]) * fe.dx(V)
                elif (par['spatialDimension'] == 2):
                    lhs += systemMatricesDomain['Kxx'][i][j] * fe.dot(fe.grad(u[j])[0], fe.grad(v[i])[0]) * fe.dx(V)
                    lhs += systemMatricesDomain['Kxy'][i][j] * fe.dot(fe.grad(u[j])[1], fe.grad(v[i])[0]) * fe.dx(V)
                    lhs += systemMatricesDomain['Kyx'][i][j] * fe.dot(fe.grad(u[j])[0], fe.grad(v[i])[1]) * fe.dx(V)
                    lhs += systemMatricesDomain['Kyy'][i][j] * fe.dot(fe.grad(u[j])[1], fe.grad(v[i])[1]) * fe.dx(V)
                else:
                    raise ValueError('Only spatial dimensions <= 2 implemented!')
                lhs += systemMatricesDomain['S'][i][j] * u[j] * v[i] * fe.dx(V)

        return lhs

    def defineBoundarySystem():
        """Define system of equations of the weak formulation of 
        the second-orderformulation of the PN equations on the boundary.
        
        The Boundary IDs are ordered like the entries in systemBoundary.
        """
        systemBoundary = sio.loadmat('%s_N_%d.mat' % (boundaryFilePrefix, N_PN), squeeze_me=True)
        outerNormalVectors3D = outerNormalVectors3DMat2List(systemBoundary['outerNormalVectors3D'])
        nBoundary = len(systemBoundary['directedFlux'])
        ub = computeBoundarySourceOddMoments(outerNormalVectors3D, par['boundarySource'], N_PN, par['rho'], par['spatialDimension'])
        
        integralsRobinBoundaryLhs = []
        integralsRobinBoundaryRhs = []
        for k in range(nBoundary):
            He = systemBoundary['He1'][k] - par['rho'][k] * systemBoundary['He2'][k]
            Ho = systemBoundary['Ho1'][k] - par['rho'][k] * systemBoundary['Ho2'][k]
            if np.isscalar(Ho):
                tmpa = np.linalg.solve(np.reshape(Ho, (1, 1)), np.reshape(He, (1,1)))
                tmpg = np.linalg.solve(np.reshape(Ho, (1, 1)), ub[k])
            else:
                tmpa = np.linalg.solve(Ho, He)
                tmpg = np.linalg.solve(Ho, ub[k])
                
            a = np.dot(systemBoundary['directedFlux'][k], tmpa)
            g = np.dot(systemBoundary['directedFlux'][k], tmpg)
            if np.isscalar(a):
                a = np.reshape(a, (1, 1))
                g = np.reshape(g, (1))

            #print('\n boundary:', k, ' \n a=', a, '\n g=', g)
            for i in range(nEvenMoments):
                for j in range(nEvenMoments):
                    integralsRobinBoundaryLhs.append(
                        a[i, j] * u[j] * v[i] * ds(k + 1))  # numbering in gmsh starts from 1
                integralsRobinBoundaryRhs.append(
                        g[i] * v[i] * ds(k + 1))

        return integralsRobinBoundaryLhs, integralsRobinBoundaryRhs
    

    # ------------------------------------------------------------------------
    nEvenMoments = getNumberOfEvenBasisFunctions(N_PN, par['spatialDimension'])
    u, v, V, uSol, uSolComponents, solAssignMod = initWeakFormulation()
    lhs = defineWeakSystemDomain()
    integralsRobinBoundaryLhs, integralsRobinBoundaryRhs = defineBoundarySystem()
    lhs += sum(integralsRobinBoundaryLhs)
    rhs = sum(integralsRobinBoundaryRhs)
    
    #-------- solve --------
    solverParameters = {"linear_solver": "gmres"}
    solverParameters = {}
    fe.solve(lhs == rhs, uSol, [], solver_parameters=solverParameters)
    solAssignMod(uSolComponents, uSol)

    return uSolComponents
