% CHECKPN2NDTESTCASE1 Check the formulas of the second-order formulation
%   for test case 1.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

clear

%%
syms mu z phi SIGMAA SIGMAS
spatialDimensionFull = 3;
spatialDimension = 1;

for N = [1, 3]
    fprintf('\n\nP%d\n---------------------------------------------\n', N)
    load(sprintf('../build/testCase1/systemBoundary_testCase1_N_%d.mat', N))
    %% get flux matrices
    [realSphericalHarmonicsSym] = getRealSphericalHarmonicsSym(N);
    nMomentsFull = getNumberOfBasisFunctions(N, spatialDimensionFull);
    [l, m] = linearIdx2DegOrder(0 : nMomentsFull - 1, spatialDimensionFull);
    idxRed = getReducedIdx(l, m, spatialDimension);
    b = cell2sym(realSphericalHarmonicsSym(idxRed));
    nMoments = getNumberOfBasisFunctions(N, spatialDimension);
    f = symfun(mu * (b * b.'), [mu, phi]);
    Tz = int(int(f, mu, -1, 1), phi, 0, 2 * pi);

    %% PN2nd reduction
    idxEven = linearIdxOfEvenBasis(N, spatialDimension) + 1 ;
    idxOdd = linearIdxOfOddBasis(N, spatialDimension) + 1;
    Tzee = Tz(idxEven, idxEven);
    Tzoo = Tz(idxOdd, idxOdd);
    Tzeo = Tz(idxEven, idxOdd);
    Tzoe = Tz(idxOdd, idxEven);
    
    C = -(SIGMAA + SIGMAS) * eye(nMoments);
    C(1,1) = C(1,1) + SIGMAS; % isotropic scattering
    Coo = C(idxOdd, idxOdd);
    invCoo = Coo\eye(length(idxOdd));
    Kzz = Tzeo * invCoo * Tzoe;
    fprintf('\ndomain: Kzz:\n')
    pretty(Kzz)
    
    fprintf('\ndomain: Cee:\n')
    pretty(C(idxEven, idxEven))
    
    %% boundary half moments left (z=0)
    idxOdd = linearIdxOfOddBasis(N, spatialDimension) + 1;
    IbLeft = sym(1 / 4) / pi;
    uL=  2 * pi * int(IbLeft * b(idxOdd), mu, 0, 1); % inward pointing, mu > 0 at left boundary
    
    %% boundary half moments right (z=1)
    uR = sym(zeros(length(idxOdd), 1));
    
    %% boundary system left (reflectivity = 0)
    f = symfun(b(idxOdd) * b(idxOdd).', [mu, phi]);
    HoL = int(int(f, mu, 0, 1), phi, 0, 2 * pi); % inward pointing, mu > 0 at left boundary
    
    f = symfun(b(idxOdd) * b(idxEven).', [mu, phi]);
    HeL = int(int(f, mu, 0, 1), phi, 0, 2 * pi); % inward pointing, mu > 0 at left boundary
    
    % n = [0;0;-1], n*Omega = -mu
    f = symfun(-mu * b(idxEven) * b(idxOdd).', [mu, phi]);
    FL = int(int(f, mu, -1, 1), phi, 0, 2 * pi);
    
    BlL = FL * (HoL\eye(length(idxOdd))) * HeL;
    BrL = FL * (HoL\eye(length(idxOdd)));
    fprintf('\nleft boundary (z=0): Bl\n')
    pretty(BlL)
    fprintf('\nleft boundary (z=0): Br\n')
    pretty(BrL)
    
    tmpl = BlL - directedFlux{1} * inv(Ho1{1}) * He1{1};
    tmpr = BrL - directedFlux{1} * inv(Ho1{1});
    if norm(tmpl(:), 'inf') > 1e-14 || norm(tmpr(:), 'inf') > 1e-14
       error('Formulae do not coincide with numeric result!') 
    end
    
    %% boundary system right (reflectivity = 0)
    f = symfun(b(idxOdd) * b(idxOdd).', [mu, phi]);
    HoR = int(int(f, mu, -1, 0), phi, 0, 2 * pi); % inward pointing, mu < 0 at right boundary
    
    f = symfun(b(idxOdd) * b(idxEven).', [mu, phi]);
    HeR = int(int(f, mu, -1, 0), phi, 0, 2 * pi); % inward pointing, mu < 0 at right boundary
    
    % n = [0;0;1], n*Omega = mu
    f = symfun(mu * b(idxEven) * b(idxOdd).', [mu, phi]);
    FR = int(int(f, mu, -1, 1), phi, 0, 2 * pi);
    
    BlR = FR * (HoL\eye(length(idxOdd))) * HeR;
    BrR = FR * (HoL\eye(length(idxOdd)));
    fprintf('\nright boundary (z=1): Bl\n')
    pretty(BlR)
    fprintf('\nright boundary (z=1): Br\n')
    pretty(BrR)
    
    tmpl = BlR - directedFlux{2} * inv(Ho1{2}) * He1{2};
    tmpr = BrR - directedFlux{2} * inv(Ho1{2});
    if norm(tmpl(:), 'inf') > 1e-14 || norm(tmpr(:), 'inf') > 1e-14
       error('Formulae do not coincide with numeric result!') 
    end
end