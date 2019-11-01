function [radiativeEnergy] = radiativeEnergyPNOrigIsotropicKernel1D(par, N, nGridPoints, varargin)
% RADIATIVEENERGYPNORIGISOTROPICKERNEL1D Evaluate the solution of the
%   original PN equations for isotropic scattering in 1D at given grid 
%   points. 
% 
%   For homogeneous absorption and scattering coefficient the solution can
%   be given analytically.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

if ~strcmp(par.kernel.name, 'isotropic')
   error('Only isotropic kernel 1/4/pi allowed.')
end
if par.spatialDimension > 1
   error('Only 1D!') 
end
%%
zGrid = linspace(par.meshZMin, par.meshZMax, nGridPoints);
rhoZMin = par.reflectivity(1);
rhoZMax = par.reflectivity(2);
sourceZMin = par.externalSource{1};
sourceZMax = par.externalSource{2};
outerNormalZMin = [0, 0, -1]';
outerNormalZMax = [0, 0, 1]';
homogeneousFlag = false;
if nargin == 4
   if strcmp(varargin{1}, 'homogeneous')
      homogeneousFlag = true; 
   end
end

%% ---- Generate PN system + boundary conditions----
% T_z * (d/dz u) = C * u
% left hand side
[~, ~, Tz] = fluxPNSphHarm(N, par.spatialDimension);
% right hand side
z = sym('z');
C = sym(zeros(N + 1, N + 1));
for l = 0 : N
    C(l + 1, l + 1) = - par.sigma_a(z);
    if l > 0
        C(l + 1, l + 1) = C(l + 1, l + 1) - par.sigma_s(z);
    end
end
bcZMin = getBoundaryCondition(outerNormalZMin, N, rhoZMin);
bcZMax = getBoundaryCondition(outerNormalZMax, N, rhoZMax);

%% ---- Diagonalize----
[V, D] = eig(Tz);

%ODE system in w:
% V^-1(V * D * V^-1) u_z = V^-1 * C * u => D * w_x = V^-1C * V * w
% => w_x = D^-1 * V^-1 * C * V * w = H * w

%% ---- Calculate half integrals at boundary ----
maxExactDegree = N + 20; % degree of source not known
[weights, mu, phi] = sphericalQuadratureHalf(outerNormalZMin, maxExactDegree);
vx = zeros(size(mu));
vy = zeros(size(mu));
vz = mu;
spatialDimension = 1;
bAtQuad = evalONB([mu; phi], N, spatialDimension);
idxOdd = linearIdxOfOddBasis(N, spatialDimension) + 1;
uZMin = sum((weights .* sourceZMin(vx, vy,  vz)) .* bAtQuad(idxOdd, :), 2);

[weights, mu, phi] = sphericalQuadratureHalf(outerNormalZMax, maxExactDegree);
vx = zeros(size(mu));
vy = zeros(size(mu));
vz = mu;
spatialDimension = 1;
bAtQuad = evalONB([mu; phi], N, spatialDimension);
idxOdd = linearIdxOfOddBasis(N, spatialDimension) + 1;
uZMax = sum((weights .* sourceZMax(vx, vy, vz)) .* bAtQuad(idxOdd, :), 2);

%% ---- Finding initial conditions of ODE ----
%Find initial conditions using BC
% wL = w(0) (-> wL = sym('w',[N+1,1]))
% d/dz w = H(z) * w ( wR = expm(H * (zMax - zMin)) * wL in homogeneous case)
% wR = w(zMax)
% uL = V * wL
% uR = V * wR = V * expm(S) * wL
% bcZMin * uL = (1 - rho) * uZMin
% bcZMax * uR = (1 - rho) * uZMax
H = D \ (V \ (C * V));
if homogeneousFlag
    % H does not depend on z, thus H(z) commutes with int_zMin^zMax H dz
    % and we get solution as expm(int_zMin^zMax H dz) * wL
    B = H * (z - par.meshZMin);
    S = double(subs(B, z, par.meshZMax));
    wL = [bcZMin * V; bcZMax * V * expm(S)] \ [(1 - rhoZMin) * uZMin; (1 - rhoZMax) * uZMax];
else
    % in general H(z) does not commute with int_zMin^zMax H dz, thus we
    % need to solve the ODE numerically
    odefun = @(x, y) double(subs(H, z, x)) * y;
    F = zeros(N + 1, N + 1, nGridPoints);
    y0 = @(i) double((1 : N + 1)'==i); 
    options = struct();
    if nargin == 4
        if isstruct(varargin{1})
            options = varargin{1};
        end
    end
    for i = 1 : N + 1
        [~, y] = ode45(odefun, zGrid, y0(i), options);
        F(:, i, :) = y';
    end
    wL = sym('w', [N + 1, 1]);
    wR = F(:, :, end) * wL;
    uL = V * wL;
    uR = V * wR;

    sol = solve([bcZMin * uL - (1 - rhoZMin) * uZMin == 0;
                 bcZMax * uR - (1 - rhoZMax) * uZMax == 0], wL);
    wL = double(struct2array(sol))';
end

%% ---- Calculating solution ----
%Right now we are only interested in u0
u0 = zeros(size(zGrid));
%unfortunately expm is not supported to be vectorized
for i = 1 : nGridPoints
    if homogeneousFlag
        u = V * (expm(double(subs(B, z, zGrid(i)))) * wL);
    else
        u = V * (F(:, :,i ) * wL);
    end
    u0(i) = u(1);
end
radiativeEnergy = u0 * sqrt(4 * pi);
end

function bc = getBoundaryCondition(outerNormal, N, rho)
    maxExactDegree = 2 * N;
    [weights, mu, phi] = sphericalQuadratureHalf(outerNormal, maxExactDegree);
    spatialDimension = 1;
    bAtQuad = evalONB([mu; phi], N, spatialDimension);
    idxOdd = linearIdxOfOddBasis(N, spatialDimension) + 1;
    H1 = integralOuterProduct(bAtQuad(idxOdd, :), bAtQuad, weights);

    [muRefl, phiRefl] = reflectedPoints(outerNormal, mu, phi);
    bAtQuadRefl = evalONB([muRefl; phiRefl], N, spatialDimension);
    H2 = integralOuterProduct(bAtQuad(idxOdd, :), bAtQuadRefl, weights);

    bc = H1 - rho * H2;
end

function [muRefl, phiRefl] = reflectedPoints(outerNormalVector3D, mu, phi)
    Omega = sphereParam2Cartesian(mu, phi);
    OmegaRefl = reflectionAtSurface(outerNormalVector3D, Omega);
    [muRefl, phiRefl] = sphereCartesian2Param(OmegaRefl);
end

