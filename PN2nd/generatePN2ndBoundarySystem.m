function [systemBoundary] = generatePN2ndBoundarySystem(outerNormalVectors3D, N, spatialDimension)
% GENERATEPN2NDBOUNDARYSYSTEM Generate the system matrices for the Marshak 
%   boundary conditions for the second-order formulation of the PN
%   equations.
%
%   We obtain a system of boundary conditions for each outer normal vector.
%
%   Compute Ho(1), Ho(2), He(1), He(2), <nOmega bb^T> here and invert it
%   later, otherwise we would need to fix the reflectivity paramter at this
%   point.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

if size(outerNormalVectors3D, 1) ~= 3
    error('Wrong input dimension: needs to be 3 x N.')
end
nVec = size(outerNormalVectors3D ,2);
idxEven = linearIdxOfEvenBasis(N, spatialDimension) + 1;
idxOdd = linearIdxOfOddBasis(N, spatialDimension) + 1;
systemBoundary.He1 = cell(nVec, 1);
systemBoundary.He2 = cell(nVec, 1);
systemBoundary.Ho1 = cell(nVec, 1);
systemBoundary.Ho2 = cell(nVec, 1);
systemBoundary.directedFlux = cell(nVec, 1);
for k = 1 : nVec
    n = outerNormalVectors3D(:, k);
    if abs(norm(n) - 1) > 1e-15
       error('Outer normal vectors should be normalized!') 
    end
    
    maxExactDegree = 2 * N;
    [weights, mu, phi] = sphericalQuadratureHalf(n, maxExactDegree);
    bAtQuad = evalONB([mu; phi], N, spatialDimension);
    systemBoundary.He1{k} = integralOuterProduct(bAtQuad(idxOdd, :), bAtQuad(idxEven, :), weights);
    systemBoundary.Ho1{k} = integralOuterProduct(bAtQuad(idxOdd, :), bAtQuad(idxOdd, :), weights);
    
    [muRefl, phiRefl] = reflectedPoints(n, mu, phi);
    bAtQuadRefl = evalONB([muRefl; phiRefl], N, spatialDimension);
    systemBoundary.He2{k} = integralOuterProduct(bAtQuad(idxOdd, :), bAtQuadRefl(idxEven, :), weights);
    systemBoundary.Ho2{k} = integralOuterProduct(bAtQuad(idxOdd, :), bAtQuadRefl(idxOdd, :), weights);
    
    maxExactDegree = 2 * N + 1;
    [weightsFull, muFull, phiFull] = sphericalQuadratureFull(maxExactDegree);
    bAtQuadFull = evalONB([muFull; phiFull], N, spatialDimension);
    OmegaFull = sphereParam2Cartesian(muFull, phiFull);
    systemBoundary.directedFlux{k} = integralOuterProduct(((n' * OmegaFull) .* bAtQuadFull(idxEven, :)), bAtQuadFull(idxOdd, :), weightsFull);
end

systemBoundary.outerNormalVectors3D = outerNormalVectors3D;

end

function [muRefl, phiRefl] = reflectedPoints(outerNormalVector3D, mu, phi)
    Omega = sphereParam2Cartesian(mu, phi);
    OmegaRefl = reflectionAtSurface(outerNormalVector3D, Omega);
    [muRefl, phiRefl] = sphereCartesian2Param(OmegaRefl);
end