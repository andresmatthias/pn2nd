% CHECKFORMULAETESTCASE1 Check formulas regarding test case 1.
%
%   rotational symmetry around z-Axis:
%   mu * d_z I + I = 0
%   I(z=0) = 1 / 4 / pi (left boundary, zero reflection)
%   I(z=1) = 0          (right boundary, zero reflection)
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%


clear
addpath(genpath('../PN2nd/'));
addpath('../testCases/testCase1')
addpath('../testCases/')
%% define candidate for the solution of the kinetic equation
syms mu z phi
I_muPos = exp(- z / mu) / 4 / pi;
I_muNeg = 0 * mu;

fprintf('\nKinetic problem\n---------------\n')
fprintf('Equation on domain (mu >=0), should be zero:\t %s\n', simplify(diff(I_muPos, z) * mu + I_muPos));
fprintf('Equation on domain (mu < 0), should be zero:\t %s\n', simplify(diff(I_muNeg, z) * mu + I_muNeg));
fprintf('Left boundary (mu >= 0), should be 1 / 4 / pi:\t %s\n', subs(I_muPos, z, 0))
fprintf('Right boundary (mu < 0), should be 0:\t\t %s\n', subs(I_muNeg, z, 1))

fprintf('I (mu>0):\n')
pretty(I_muPos)
fprintf('I (mu<0):\n')
pretty(I_muNeg)

%% PN analytic
% Tz * u = C * u
% HL * u = uL, HR * u = uR (left (z=0), right (z=1); zero reflectivity)

fprintf('\n')
spatialDimensionFull = 3;
spatialDimension = 1;
for N = [1, 3]
    fprintf('\n\nP%d\n---------------------------------------------\n', N)
    %% flux matrix
    [~, ~, TzNum] = fluxPNSphHarm(N, spatialDimension);

    [realSphericalHarmonicsSym] = getRealSphericalHarmonicsSym(N);
    nMomentsFull = getNumberOfBasisFunctions(N, spatialDimensionFull);
    [l, m] = linearIdx2DegOrder(0 : nMomentsFull - 1, spatialDimensionFull);
    idxRed = getReducedIdx(l, m, spatialDimension);

    b = cell2sym(realSphericalHarmonicsSym(idxRed));
    nMoments = getNumberOfBasisFunctions(N, spatialDimension);
    f = symfun(mu * (b * b.'), [mu, phi]);
    Tz = int(int(f, mu, -1, 1), phi, 0, 2 * pi);
    
    tmp = TzNum - vpa(Tz);
    fprintf('Difference symbolic / numeric Tz:\nshould be close to zero:\t\t\t%1.2e \n', norm(tmp(:), 'inf'))
    fprintf('\nTz:\n\n')
    disp(Tz)
    %% boundary half moments left (z=0)
    idxOdd = linearIdxOfOddBasis(N, spatialDimension) + 1;
    IbLeft = sym(1 / 4) / pi;
    uL =  2 * pi * int(IbLeft * b(idxOdd), mu, 0, 1); % inward pointing, mu > 0 at left boundary
    
    %% boundary half moments right (z=1)
    uR = sym(zeros(length(idxOdd), 1));
    
    %% boundary system left (reflectivity = 0)
    % HL * u = uL
    f = symfun(b(idxOdd) * b.', [mu, phi]);
    HL = int(int(f, mu, 0, 1), phi, 0, 2 * pi); % inward pointing, mu > 0 at left boundary
    fprintf('\nHL:\n\n')
    pretty(HL)
    fprintf('\nuL:\n\n')
    pretty(uL)
    %% boundary system right (reflectivity = 0)
    % HR * u = uR
    idxOdd = linearIdxOfOddBasis(N, spatialDimension) + 1;
    f = symfun(b(idxOdd) * b.', [mu, phi]);
    HR = int(int(f, mu, -1, 0), phi, 0, 2 * pi); % inward pointing, mu > 0 at left boundary
    fprintf('\nHR:\n\n')
    pretty(HR)
    fprintf('\nuR:\n\n')
    pretty(uR)
    %% collision operator (zero scattering)
    C = -eye(nMoments);
     fprintf('\nC:\n\n')
    disp(C)
    %% diagonalize
    [V, D] = eig(Tz);
    A = D\(V * C * (V\eye(nMoments)));

    %% solution
    if N == 1
        % u = V * w
        % HL * u = uL; HR * u = uR
        % solve for initial condition of w
        wL = [HL * V; HR * V * expm(A)] \ [uL; uR];
        wL = simplify(wL);

        u = V * expm(A * z) * wL;
        u = simplify(u);

        %% validate with equation
        fprintf('Equation on domain, should be zero:\t\t %s\n', simplify(Tz * diff(u, z) - C * u));

        %% validate with test case 1 script
        warning('off', 'MATLAB:mpath:nameNonexistentOrNotADirectory')
        par = loadTestCase1();
        warning('off', 'MATLAB:mpath:nameNonexistentOrNotADirectory')
        nGridPoints = 50;
        [radiativeEnergy] = radiativeEnergyPNOrigIsotropicKernel1D(par, N, nGridPoints, 'homogeneous');
        radiativeEnergy2 = sqrt(4 * pi) * subs(u(1), z, linspace(0, 1, nGridPoints));
        fprintf('Difference between two solutions\nshould be close to zero:\t\t\t%e\n', norm(radiativeEnergy2 - radiativeEnergy, 'inf'))

        fprintf('\nu:\n\n')
        pretty(u)
    end
end
