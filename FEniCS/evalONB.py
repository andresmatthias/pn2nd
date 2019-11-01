"""Evaluate orthonormal basis on sphere (real spherical harmonics).

For details on how to use the defined functions, see
    ../unitTests.py

%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
"""

import numpy as np
from scipy.special import factorial

def evalONB(pointsOnSphereMuPhi, N, spatialDimension):
    """Evaluate ONB on sphere. Depending on spatialDimension we can reduce
    the number of basis functions due to symmetry assumptions.
    """
    basisAtPoints = evalRealSphericalHarmonics(pointsOnSphereMuPhi, N)
    fullSpatialDimension = 3
    nMomentsFull = getNumberOfBasisFunctions(N, fullSpatialDimension)
    linIdx = np.arange(0, nMomentsFull)
    l, m = linearIdx2DegOrder(linIdx, fullSpatialDimension)
    idxReduced = getReducedIdx(l, m, spatialDimension) 
    basisAtPoints = basisAtPoints[idxReduced, :] 
    
    return basisAtPoints
    
    
def evalRealSphericalHarmonics(pointsOnSphereMuPhi, lmax):
    """Evaluate real spherical harmonics up to degree lmax, based on
    
    Evaluation of rotation matrices in the basis of real spherical
    harmonics; Miguel A. Blanco, M. Florez, M. Bermejo; 1997; Journal
    of Molecular Structure: THEOCHEM; Volume 419; 
    doi: https://doi.org/10.1016/S0166-1280(97)00185-1
    """
    P = evalLegendreAssociatedPolynomials(lmax, pointsOnSphereMuPhi)
    C = [np.cos(pointsOnSphereMuPhi[1, :]) for m in range(0, lmax)]
    S = [np.sin(pointsOnSphereMuPhi[1, :]) for m in range(0, lmax)]

    for m in range(1, lmax):
        C[m] = C[0] * C[m - 1] - S[0] * S[m - 1]
        S[m] = S[0] * C[m - 1] + C[0] * S[m - 1]

    spatialDimension = 3
    nMoments =  getNumberOfBasisFunctions(lmax, spatialDimension)
    realSphericalHarmonicsAtPoints = np.zeros((nMoments, pointsOnSphereMuPhi.shape[1]))
    N = lambda l, m: np.sqrt((2 * l + 1) / 4 / np.pi * factorial(l - m) / factorial(l + m))
    linIdx = np.arange(0, nMoments)
    l = np.floor(np.sqrt(linIdx) + 1e-15).astype(int)
    m = linIdx - l ** 2 - l
    for i in range(0, len(linIdx)):
        if m[i] < 0:
            realSphericalHarmonicsAtPoints[linIdx[i], :] = np.sqrt(2) * S[abs(m[i]) - 1] * N(l[i], abs(m[i])) * P[l[i]][abs(m[i])]
        if m[i] == 0:
            realSphericalHarmonicsAtPoints[linIdx[i], :] = N(l[i], m[i]) * P[l[i]][m[i]]
        if m[i] > 0:
            realSphericalHarmonicsAtPoints[linIdx[i], :] = np.sqrt(2) * C[m[i] - 1] * N(l[i], m[i]) * P[l[i]][m[i]]

    return realSphericalHarmonicsAtPoints


def evalLegendreAssociatedPolynomials(lmax, pointsOnSphereMuPhi):
    """Evaluate Legendre associated polynomials up to degree lmax, based on
    
    Evaluation of rotation matrices in the basis of real spherical
    harmonics; Miguel A. Blanco, M. Florez, M. Bermejo; 1997; Journal
    of Molecular Structure: THEOCHEM; Volume 419; 
    doi: https://doi.org/10.1016/S0166-1280(97)00185-1
        
    Note: Implementation without Condon-Shortley phase! 
    """
    mu = pointsOnSphereMuPhi[0, :]
    P = [[None for m in range(lmax + 1)] for l in range(lmax + 1)]
    l = 0
    P[l][l] = np.ones(mu.shape)
    for l in range(1, lmax + 1):
        P[l][l] = (2 * l - 1) * np.sqrt(1 - mu**2) * P[l - 1][l - 1]
        P[l][l - 1] = (2 * (l - 1) + 1) * mu * P[l - 1][l - 1]
    for l in range(2, lmax + 1):
        for m in range(0, l - 1):
            P[l][m] = ( (2 * l - 1) * mu * P[l - 1][m] - (l + m - 1) * P[l - 2][m] ) / (l - m)

    return P


def linearIdx2DegOrder(linIdx, spatialDimension):
    """ 3D (-l<=m<=l)     2D(l+m even)  1D(m=0)
        n: (l, m)         n: (l, m)     n: (l, m)
        0: (0, 0)         0: (0, 0)     0: (0, 0)
        1: (1, -1)        1: (1, -1)    
        2: (1, 0)                       1: (1, 0)
        3: (1, 1)         2: (1, 1)
        4: (2, -2)        3: (2, -2)
        5: (2, -1)        
        6: (2, 0)         4: (2, 0)     2: (2, 0)
        7: ...
    """
    
    if spatialDimension == 1:
        l = linIdx
        m = np.zeros(l.shape, dtype='int')
    elif spatialDimension == 2:
        l = np.ceil(-3 / 2 + np.sqrt(9 / 4 + 2 * linIdx) - 1e-15).astype('int')
        m = (2 * (linIdx - l * (l + 1) / 2) - l).astype('int')
    elif spatialDimension == 3:
        l = np.floor(np.sqrt(linIdx) + 1e-15).astype('int')
        m = (linIdx - l**2 - l).astype('int')
    else:
        raise ValueError('Invalid spatial dimension')
    
    return l, m


def getNumberOfBasisFunctions(N, spatialDimension):
    if spatialDimension == 1:
        return N + 1
    elif spatialDimension == 2:
        return N ** 2 / 2 + 3 / 2 * N + 1
    elif spatialDimension == 3:
        return N ** 2 + 2 * N + 1
    else:
        raise ValueError('Invalid spatial dimension')


def linearIdxOddBasis(N, spatialDimension):
    """Compute linear indices of odd basis function, (those with odd degree)."""
    nMoments = getNumberOfBasisFunctions(N, spatialDimension)
    linIdx = np.arange(0, nMoments)
    l, _ = linearIdx2DegOrder(linIdx, spatialDimension)
    idxOdd = np.where(np.mod(l, 2) == 1)
    return idxOdd


def getReducedIdx(l, m, spatialDimension):
    """Return indices of those basis functions included after reduction of 
    dimension due to symmetry assumptions.
    """
    if spatialDimension == 1:
        idxRed = (m == 0)
    elif spatialDimension == 2:
        idxRed = (np.mod(l + m, 2) == 0)
    elif spatialDimension == 3:
        idxRed = np.ones(l.shape, dtype='bool')
    
    return idxRed
