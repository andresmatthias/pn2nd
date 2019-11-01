"""Compute the odd half moments at the boundary (integral over inward
pointing directions; only take basis functions (real spherical harmonics) with
odd degree).

For details on how to use the defined functions, see
    ../unitTests.py

%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
"""

import numpy as np
import scipy.io as sio
import evalONB as basis
import os


def computeBoundarySourceOddMoments(outerNormalVectors, boundarySource, lmax, rho, spatialDimension):
    """S(vx, vy, vz), can be evaluated vectorized"""
    ubOdd = []
    idxOdd = basis.linearIdxOddBasis(lmax, spatialDimension)
    for k in range(len(outerNormalVectors)):
        mu, phi, weights = getSphericalQuadratureHalf(outerNormalVectors[k])
        basisAtQuad = basis.evalONB(np.vstack((mu, phi)), lmax, spatialDimension)

        Omega = sphereParam2Cartesian(mu, phi)
        boundarySourceAtQuad = boundarySource[k](Omega[0, :], Omega[1, :], Omega[2, :])
        boundarySourceAtQuad = np.expand_dims(boundarySourceAtQuad, axis=0)
        boundarySourceAtQuad = np.tile(boundarySourceAtQuad, (basisAtQuad.shape[0], 1))
        weights = np.tile(weights, (basisAtQuad.shape[0], 1))
        allMoments = np.sum(weights * (1 - rho[k]) * boundarySourceAtQuad * basisAtQuad , axis=1)
        ubOdd.append(allMoments[idxOdd])

    return ubOdd


def getSphericalQuadratureHalf(outerNormalVector):
    scriptDir = os.path.dirname(__file__)
    relPath = '/sphericalQuadratureUpperHalf_ExactDegree_50.mat'
    QUpper = sio.loadmat(scriptDir + relPath)  
    Omega0 = sphereParam2Cartesian(QUpper['mu'], QUpper['phi'])
    R = getGeomRotationa2b(np.array([0, 0, -1]), outerNormalVector)
    Omega = np.dot(R, Omega0)
    mu, phi = sphereCartesian2Param(Omega)
    
    return mu, phi, QUpper['weights']


def getGeomRotationa2b(a, b):
    """Geometrical rotation of vector a to b."""
    if (np.abs(np.linalg.norm(a) - 1) + np.abs(np.linalg.norm(b) - 1)) > 1e-15:
        raise ValueError('Vectors have to be normalized!')
        
    if np.linalg.norm(a + b) < 1e-15:
    # formula is not applicable in this case but integrals can be
    # computed analogously
        return np.diag([1,1,-1])
    else:
        v = np.cross(a, b)
        c = np.dot(a, b)
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]],[ -v[1], v[0], 0]])
        return (np.eye(3) + vx + np.dot(vx, vx) / (1 + c))


def sphereParam2Cartesian(mu, phi):
    """Get cartesian representation of point on sphere."""
    Omega = np.vstack((np.sqrt(1 - mu**2) * np.cos(phi),
                       np.sqrt(1 - mu**2) * np.sin(phi),
                       mu))
    
    return Omega


def sphereCartesian2Param(pointsOnSphere):
    """Get parametrization in z-coordinate and polar angle of point on sphere"""
    mu = pointsOnSphere[2, :]
    phi = np.arctan2(pointsOnSphere[1, :], pointsOnSphere[0, :])
    
    return mu, phi
