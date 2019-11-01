% FLUXMATRIXTEST Perform unit tests regarding properties of the flux 
%   matrices of the PN equations, depending on symmetry assumptions
%   (spatialDimension)).
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

addpath(genpath('./quadratureRules/'))
addpath('./basisFunctions/')


%% explicit and numerical PN fluxes coincide 
for spatialDimension = [1, 2, 3]
    N = 10;
    maxExactDegree = 2 * N + 1; % basis_i * basis_j * Omega_?, deg(basis) <= N
    [weights, mu, phi] = sphericalQuadratureFull(maxExactDegree);
    b = evalONB([mu; phi], N, spatialDimension);

    nMoments = size(b, 1);
 
    Omega = sphereParam2Cartesian(mu, phi);
    Omegax =  repmat(Omega(1, :), nMoments, 1);
    Omegay =  repmat(Omega(2, :), nMoments, 1);
    Omegaz =  repmat(Omega(3, :), nMoments, 1);
    TxNum = integralOuterProduct(Omegax .* b, b, weights);  % is of degree 2*N + 1
    TyNum = integralOuterProduct(Omegay .* b, b, weights);  % is of degree 2*N + 1
    TzNum = integralOuterProduct(Omegaz .* b, b, weights);  % is of degree 2*N + 1

    [Tx, Ty, Tz] = fluxPNSphHarm(N, spatialDimension);
    assert(norm(Tx - TxNum, 'inf') < 1e-12);
    assert(norm(Ty - TyNum, 'inf') < 1e-12);
    assert(norm(Tz - TzNum, 'inf') < 1e-12);
end

