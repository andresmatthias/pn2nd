% PN2ndREDUCTIONTEST Perform unit tests regarding the second-order
% formulation of the PN equations.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

%% projection of scattering operator
addpath('./basisFunctions/')
addpath('./quadratureRules/')
for spatialDimension = [1, 2, 3]
    N = 10;
    % this is not a scattering kernel, but we use it to validate computation of
    % the double integral
    kernelFun = @(vx, vy, vz, wx, wy, wz) ones(size(vx .* wx)); 
    CScatter = projectScatteringOperator(N, kernelFun, spatialDimension);
    CScatterTrue = zeros(size(CScatter));
    CScatterTrue(1,1) = 4 * pi;
    assert(norm(CScatter - CScatterTrue, 'inf') < 1e-12)
end