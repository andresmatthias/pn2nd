"""Wrappers for 1D and 2D solution of the second-order formulation of the PN 
equations.

%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%
"""

import os
import fenics as fe
import numpy as np
import scipy.io as sio
from solvePN2nd import solvePN2nd


def wrapper2D(par, meshName, N_PN, boundaryFilePrefix, solutionFolder):
    #-------- load mesh --------
    mesh = fe.Mesh(meshName + '.xml')
    boundaryMarkers = fe.MeshFunction('size_t', mesh, meshName + '_facet_region.xml')
    ds = fe.Measure('ds', domain=mesh, subdomain_data=boundaryMarkers)
    
    u = solvePN2nd(par, N_PN, mesh, boundaryFilePrefix, ds)
    #fe.plot(u[0])

    # first moment:     u0 = int_{4pi} I * b0 dOmega, b0 = sqrt(1 / 4 / pi)
    # radiative energy: phi = int_{4pi} I dOmega 
    radiativeEnergy = np.sqrt(4 * np.pi) * u[0].compute_vertex_values()
    createFolder(solutionFolder)
    sio.savemat('%sradiativeEnergy_P%d2nd_%s.mat'%(solutionFolder, N_PN, meshName.split('/')[-1]),
                    {'radiativeEnergy': radiativeEnergy, 'points': mesh.coordinates(),
                     'connectivityList': mesh.cells()})
    
    return u



def wrapper1D(par, nGridPoints, N_PN, boundaryFilePrefix, solutionFolder): 
    #-------- load mesh --------
    mesh = fe.IntervalMesh(nGridPoints - 1, par['domain']['zMin'], par['domain']['zMax'])
    ds = defineBoundaryMarkers1D(mesh, par['domain']['zMin'], par['domain']['zMax'])
    
    u = solvePN2nd(par, N_PN, mesh, boundaryFilePrefix, ds)
    #fe.plot(u[0])
    
    # first moment:     u0 = int_{4pi} I * b0 dOmega, b0 = sqrt(1 / 4 / pi)
    # radiative energy: phi = int_{4pi} I dOmega 
    radiativeEnergy = np.sqrt(4 * np.pi) * u[0].compute_vertex_values()
    createFolder(solutionFolder)
    sio.savemat('%sradiativeEnergy_P%d2nd_n_%d.mat'%(solutionFolder, N_PN, nGridPoints),
                    {'radiativeEnergy': radiativeEnergy, 'points': mesh.coordinates(),
                     'connectivityList': mesh.cells()})
    
    return u


def defineBoundaryMarkers1D(mesh, zMin, zMax):
    """For implementation details, see
    https://fenicsproject.org/olddocs/dolfin/1.3.0/python/demo/documented/subdomains-poisson/python/documentation.html
    https://fenicsproject.discourse.group/t/natural-boundary-conditions-without-facetfunction/228/3
    """
    
    boundaries = fe.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    class Left(fe.SubDomain):
        def inside(self, x, on_boundary):
            return fe.near(x[0], zMin)

    class Right(fe.SubDomain):
        def inside(self, x, on_boundary):
            return fe.near(x[0], zMax)
    
    left = Left()
    right = Right()
    left.mark(boundaries, 1)
    right.mark(boundaries, 2)
    ds = fe.Measure("ds", domain=mesh, subdomain_data=boundaries)
    
    return ds


def createFolder(folder):
    if not os.path.exists(folder):
        os.mkdir(folder)
