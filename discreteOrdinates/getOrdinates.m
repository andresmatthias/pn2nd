function ordinates = getOrdinates(maxExactDegree, dimension, kernelName)
% GETORDINATES Get discrete ordinates from Gaussian like quadrature rules 
%   on unit sphere. Depending on symmetry assumptions we can restrict 
%   ourselves to a subset of those discrete  ordinates.
% 
%   For 1D isotropic: only take ordinates with one specific phi. 
%   Our choice of quadrature nodes excludes mu = 0 in the 1D case, which
%   would cause a zero column in the transport term.
% 
%   In 2D isotropic: ony take ordinates with mu >= 0.
%   Our choice of quadrature nodes excludes mu = 1 in the 2D case, which 
%   would cause a zero column in the transport term.
% 
%   For a demonstration on how to use this function,
%   see also MISCDISCRETEORDINATESTEST
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

if strcmp(kernelName, 'isotropic') && dimension == 1
    % to avoid mu = 0, concatenate quadrature on upper and lower half
    [weights, mu, phi] = sphericalQuadratureUpperHalf(maxExactDegree);
    weights = [weights, weights];
    mu = [mu, -mu];
    phi = [phi, phi];
    
    tmp = abs(phi' - unique(phi));
    tmp = unique(tmp(:));
    minDiffPhi = tmp(2);
    
    idx = abs(phi - phi(1)) < minDiffPhi / 10;
    
    if mod(length(phi), sum(idx)) ~= 0
       error('Selection of phi went wrong') 
    end
    weights = weights(idx) * length(unique(phi));
    mu = mu(idx);
    phi = 0 * phi(idx);
    
elseif strcmp(kernelName, 'isotropic') && dimension == 2
    [weights, mu, phi] = sphericalQuadratureUpperHalf(maxExactDegree);
    weights = 2 * weights;
    
else
% for anisotropic kernel take full sphere, e.g., for Henyey-Greenstein.
    [weights, mu, phi] = sphericalQuadratureFull(maxExactDegree);
end

ordinates.weights = weights;
ordinates.mu = mu;
ordinates.phi = phi;
ordinates.Cartesian = sphereParam2Cartesian(ordinates.mu, ordinates.phi);
end
