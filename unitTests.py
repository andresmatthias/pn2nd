"""Execute unit tests in the context of solving the second-order formulation
of the PN equations.

%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%
"""


import unittest
import sympy as sp
import numpy as np
from scipy.special import factorial
import scipy.io as sio
import scipy.integrate as integrate
import sys
PN2ndFEniCSFolder = './FEniCS'
sys.path.append(PN2ndFEniCSFolder)
import evalONB as basis
import boundarySourceOddMoments as boundaryMoments


class TestBasisFunctions(unittest.TestCase):
    """Real spherical harmonics as basis functions on sphere."""

    def test_LegendreAssociatedPolynomials(self):
        """Compare stable implementation of associated Legendre polynomials with
        Rodrigues' formula, without Condon-Shortley phase!"""
        mu = sp.Symbol('mu')
        muGrid = np.linspace(-1, 1, 100)
        lmax = 5
        P = basis.evalLegendreAssociatedPolynomials(lmax, np.reshape(muGrid, (1, len(muGrid))))
        cnt = 0
        for l in range(0, lmax + 1):
            for m in range(0, l + 1):
                cnt += 1
                PExact = 1 / (2 ** l * factorial(l)) * (1 - mu ** 2) ** (m / 2) * sp.diff((mu ** 2 - 1) ** l, mu, l + m)
                PExactEval = [PExact.subs(mu, imu) for imu in muGrid]
                self.assertLess(np.linalg.norm(np.array(PExactEval, dtype='float') - P[l][m]), 1e-11)

    def test_sumQuadratureWeightsFull(self):
        Q = sio.loadmat(PN2ndFEniCSFolder + '/sphericalQuadratureFull_ExactDegree_50.mat')
        self.assertLess(np.abs(np.sum(Q['weights']) - 4 * np.pi), 1e-13)

    def test_ONB(self):
        """Real spherical harmonics are orthonormal basis."""
        lmax = 9
        spatialDimension = 3
        Q = sio.loadmat(PN2ndFEniCSFolder + '/sphericalQuadratureFull_ExactDegree_50.mat')
        S = basis.evalONB(np.vstack((Q['mu'], Q['phi'])), lmax, spatialDimension)
        SAux = np.expand_dims(S, axis=0)
        SAux = np.tile(SAux, (SAux.shape[1], 1, 1))
        SAuxTranspose = np.expand_dims(S, axis=1)
        SAuxTranspose = np.tile(SAuxTranspose, (1, SAux.shape[1], 1))
        wAux = np.expand_dims(Q['weights'], axis=0)
        wAux = np.tile(wAux, (SAux.shape[0], SAux.shape[1], 1))
        integral = np.sum(SAux * SAuxTranspose * wAux, axis=2)
        self.assertLess(np.linalg.norm(integral - np.eye(integral.shape[0])), 1e-13)
    
    def test_linearIdx2DegOrder(self):
        """Conversion from linear index to (l, m)-tuple index."""
        linIdx = np.arange(0, 101)
        for spatialDimension in [1, 2, 3]:
            ref = sio.loadmat(PN2ndFEniCSFolder + '/unitTestReferenceFiles/linIdx1_100_to_DegOrder_d%d.mat' % spatialDimension)
            l, m = basis.linearIdx2DegOrder(linIdx, spatialDimension)
            self.assertEqual(np.linalg.norm(l - ref['l']), 0)
            self.assertEqual(np.linalg.norm(m - ref['m']), 0)
            
    def test_linearIdxOddBasis(self):
        N = 10
        for spatialDimension in [1, 2, 3]:
            ref = sio.loadmat(PN2ndFEniCSFolder + '/unitTestReferenceFiles/linIdxOdd_N10_d%d.mat' % spatialDimension)
            idxOdd = basis.linearIdxOddBasis(N, spatialDimension)
            self.assertEqual(np.linalg.norm(idxOdd - ref['idxOdd']), 0)

    def test_reducedIndices(self):
        linIdx = np.arange(0, 101)
        spatialDimension = 3
        l, m = basis.linearIdx2DegOrder(linIdx, spatialDimension)
        for spatialDimension in [1, 2, 3]:
            ref = sio.loadmat(PN2ndFEniCSFolder + '/unitTestReferenceFiles/reducedIdx_1_100_d%d.mat' % spatialDimension)
            redIdx = basis.getReducedIdx(l, m, spatialDimension)
            self.assertEqual(np.sum(ref['redIdx'] - redIdx.astype('int')), 0)


def getRealSphericalHarmonicsFun():
    """Spherical harmonics for (l,m) = (0,0), (1, -1), ..., (2, -1), (2, 2)"""
    S = [lambda mu, phi: np.sqrt(1 / 4 / np.pi)]
    S.append(lambda mu, phi: 0.48860251190291994921581851419453 * np.sin(phi) * (1.0 - 1.0 * mu ** 2) ** (1/2))
    S.append(lambda mu, phi: 0.48860251190291992262615394793102 * mu)
    S.append(lambda mu, phi: 0.48860251190291994921581851419453 * np.cos(phi) * (1.0 - 1.0 * mu ** 2) ** (1/2))
    S.append(lambda mu, phi: -0.36418281019735973495343096974567 * np.cos(phi) * np.sin(phi) * (3.0 * mu **2 - 3.0))
    S.append(lambda mu, phi: 1.092548430592079204860292909237 * mu * np.sin(phi) * (1.0 - 1.0 * mu**2) ** (1/2))
    S.append(lambda mu, phi: 0.94617469575756013577816361248551 * mu**2 - 0.3153915652525200452593878708285)
    S.append(lambda mu, phi: 1.092548430592079204860292909237 * mu * np.cos(phi) * (1.0 - 1.0 * mu**2) ** (1/2))
    S.append(lambda mu, phi: -0.18209140509867986747671548487284 * (3.0 * mu**2 - 3.0)*(np.cos(phi) ** 2 - 1.0 * np.sin(phi)**2))

    return S


class TestBoundarySourceMoments(unittest.TestCase):
    """Compute integral over inward-pointing directions of source at boundary"""
    def test_rotationMatrix(self):
        n = 50
        for k in range(0, n):
            a = np.random.randn(3)
            b = np.random.randn(3)
            a = a / np.linalg.norm(a)
            b = b / np.linalg.norm(b)
            R = boundaryMoments.getGeomRotationa2b(a, b)
            self.assertLess(np.linalg.norm(np.dot(R, a) - b), 1e-13)

    def test_parametrizationOfSphere(self):
        n = 50
        mu = np.reshape(np.linspace(-1, 1, n), (1, n))
        phi = np.reshape(np.linspace(-np.pi, np.pi, n), (1, n))
        mu1, phi1 = boundaryMoments.sphereCartesian2Param(boundaryMoments.sphereParam2Cartesian(mu, phi))
        self.assertLess(np.linalg.norm(mu1 - mu), 1e-14)
        self.assertLess(np.linalg.norm(phi1 - phi), 1e-14)

    def test_parametrizationOfSphereInv(self):
        n = 50
        v = np.random.rand(3, n)
        v = v / np.linalg.norm(v, axis=0)
        self.assertLess(np.sum(np.linalg.norm(v, axis=0) - 1), 1e-14)
        mu, phi = boundaryMoments.sphereCartesian2Param(v)
        v1 = boundaryMoments.sphereParam2Cartesian(mu, phi)
        self.assertLess(np.linalg.norm(v1 - v), 1e-14)

    def test_sumQuadratureWeightsHalf(self):
        Q = sio.loadmat(PN2ndFEniCSFolder + '/sphericalQuadratureUpperHalf_ExactDegree_50.mat')
        self.assertLess(np.abs(np.sum(Q['weights']) - 2 * np.pi), 1e-13)

    def test_quadratureHalfZgeq0(self):
        """Half space integration, z > 0."""
        S = getRealSphericalHarmonicsFun()
        for k in range(0, len(S)):
            intRef = integrate.dblquad(S[k],  0, 2 * np.pi, lambda mu: 0, lambda mu: 1.0)
            mu, phi, weights = boundaryMoments.getSphericalQuadratureHalf(np.array([0, 0, -1]))
            lmax = 2
            pointsOnSphereMuPhi = np.vstack((mu, phi))
            SNum = basis.evalRealSphericalHarmonics(pointsOnSphereMuPhi, lmax)
            intNum = np.sum(SNum[k] * weights)

            self.assertLess(np.abs(intNum - intRef[0]), 1e-14)

    def test_quadratureHalfXgeq0(self):
        """Half space integration, x > 0."""
        S = getRealSphericalHarmonicsFun()
        for k in range(0, len(S)):
            intRef = integrate.dblquad(S[k], -np.pi / 2, np.pi / 2, lambda mu: -1.0, lambda mu: 1.0)
            mu, phi, weights = boundaryMoments.getSphericalQuadratureHalf(np.array([-1, 0, 0]))
            lmax = 2
            pointsOnSphereMuPhi = np.vstack((mu, phi))
            SNum = basis.evalRealSphericalHarmonics(pointsOnSphereMuPhi, lmax)
            intNum = np.sum(SNum[k] * weights)

            self.assertLess(np.abs(intNum - intRef[0]), 1e-14)


if __name__ == '__main__':
    unittest.main()