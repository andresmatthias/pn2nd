function [systemMatricesDomain] = generatePN2ndWeakSystemOnDomain(N, kernel, spatialDimension)
% GENERATEPN2NDWEAKSYSTEMONDOMAIN Generate the second-order formulation of
%   the PN equations for the interior of the spatial domain.
%
%   The kernel function needs to allow vectorized evaluation, e.g., 
%   k = @(vx, vy, vz, wx, wy, wz) 1 / 4 / pi * ones(size(vx .* wx))
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

% addpath('./quadratureRules/')
tolCeoCoe = 1e-12; % stop, if drift term occurs
if mod(N, 2) ~= 1
   error('Only odd model orders allowed for second-order formulation PN!') 
end

%% compute scattering operator
idxEven = linearIdxOfEvenBasis(N, spatialDimension) + 1;
idxOdd = linearIdxOfOddBasis(N, spatialDimension) + 1;
nMoments = getNumberOfBasisFunctions(N, spatialDimension);
sigma_a = str2sym('sigma_a');
sigma_s = str2sym('sigma_s');
sigma_t = str2sym('sigma_t');
sigmaAux = str2sym('sigmaAux');

CScatterAux = projectScatteringOperator(N, kernel.fun, spatialDimension);
CAttenuationAux = eye(nMoments);
if (norm(CScatterAux(idxEven, idxOdd), 'inf') > tolCeoCoe) ||...
        (norm(CScatterAux(idxOdd, idxEven), 'inf') > tolCeoCoe) 
    error('C_eo or C_oe contain non-zero entries; this would lead to drift term; stop here, as no numerical solution is provided')    
end
C = -(sigma_a + sigma_s) * CAttenuationAux + sigma_s * CScatterAux;
invCoo = computeInverseCoo();

%% compute flux matrices
[Tx, Ty, Tz] = fluxPNSphHarm(N, spatialDimension);

%% generate system, without drift term
%  (Kxx*dx_uEven + Kxy*dy_uEven + Kxz*dz_uEven) * dx_phi 
% +(Kyx*dx_uEven + Kyy*dy_uEven + Kyz*dz_uEven) * dy_phi
% +(Kzx*dx_uEven + Kzy*dy_uEven + Kzz*dz_uEven) * dz_phi
% + S*uEven = 0
switch spatialDimension
    case 1
        systemMatricesDomain.Kzz = Tz(idxEven, idxOdd) * (invCoo * Tz(idxOdd, idxEven));
    case 2
        systemMatricesDomain.Kxx = Tx(idxEven, idxOdd) * (invCoo * Tx(idxOdd, idxEven));
        systemMatricesDomain.Kxy = Tx(idxEven, idxOdd) * (invCoo * Ty(idxOdd, idxEven));

        systemMatricesDomain.Kyx = Ty(idxEven, idxOdd) * (invCoo * Tx(idxOdd, idxEven));
        systemMatricesDomain.Kyy = Ty(idxEven, idxOdd) * (invCoo * Ty(idxOdd, idxEven));
    case 3
        systemMatricesDomain.Kxx = Tx(idxEven, idxOdd) * (invCoo * Tx(idxOdd, idxEven));
        systemMatricesDomain.Kxy = Tx(idxEven, idxOdd) * (invCoo * Ty(idxOdd, idxEven));
        systemMatricesDomain.Kxz = Tx(idxEven, idxOdd) * (invCoo * Tz(idxOdd, idxEven));

        systemMatricesDomain.Kyx = Ty(idxEven, idxOdd) * (invCoo * Tx(idxOdd, idxEven));
        systemMatricesDomain.Kyy = Ty(idxEven, idxOdd) * (invCoo * Ty(idxOdd, idxEven));
        systemMatricesDomain.Kyz = Ty(idxEven, idxOdd) * (invCoo * Tz(idxOdd, idxEven));

        systemMatricesDomain.Kzx = Tz(idxEven, idxOdd) * (invCoo * Tx(idxOdd, idxEven));
        systemMatricesDomain.Kzy = Tz(idxEven, idxOdd) * (invCoo * Ty(idxOdd, idxEven));
        systemMatricesDomain.Kzz = Tz(idxEven, idxOdd) * (invCoo * Tz(idxOdd, idxEven));   
end
systemMatricesDomain.S = C(idxEven, idxEven);
fn = fieldnames(systemMatricesDomain);
for k = 1:length(fn)
    if strcmp(kernel.name, 'isotropic')
        systemMatricesDomain.(fn{k}) = simplify(systemMatricesDomain.(fn{k}));
    else
        fprintf('Clean up symbolic matrix %s\n', fn{k});
        systemMatricesDomain.(fn{k}) = cleanSymbolicMatrix(systemMatricesDomain.(fn{k}));
    end
end


function invCoo = computeInverseCoo()
    CooAux = CScatterAux(idxOdd, idxOdd) - sigmaAux * CAttenuationAux(idxOdd, idxOdd);
    invCooAux = choleskySolve(CooAux, eye(size(CooAux)));
    invCooAux = subs(invCooAux, sigmaAux, sigma_t / sigma_s);
    invCoo = simplify(invCooAux / sigma_s);
end

end