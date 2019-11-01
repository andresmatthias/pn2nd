% SPHERICALQUADRATURETEST Perform unit tests regarding quadrature rules for 
%   functions defined on the unit sphere (like the real spherical   
%   harmonics). 
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

%% weights of spherical quadrature sum up to 4*pi full / half sphere
for m = 2:10:102
    [weights] = sphericalQuadratureFull(m);
    assert(abs(sum(weights) - 4 * pi) < 1e-12);
    [weights] = sphericalQuadratureUpperHalf(m);
    assert(abs(sum(weights) - 2 * pi) < 1e-12);
end

%% trigonometric gauss
for k = 2 : 20
   % trigonometric polynomial of degree k
   f = @(phi) cos(k * phi);
   g = @(phi) sin(k * phi);
   tmp = rand(1);
   % quadrature rule should be exact for degree k
   tw = trigauss(k, tmp, tmp * 2 * pi);
   assert(abs(integral(f, tmp, tmp * 2 * pi) - sum(f(tw(:, 1)) .* tw(:, 2))) < 1e-13);
end

%% quadrature on full sphere
addpath('../basisFunctions/')
mu = sym('mu');
phi = sym('phi'); 
Nmax = 5;
for N = 2 : Nmax  % no rule available for N = 1
    realSphericalHarmonics = getRealSphericalHarmonicsSym(N);
    [weights, muQuad, phiQuad] = sphericalQuadratureFull(N);
    for k = length(realSphericalHarmonics) - 2 * N : length(realSphericalHarmonics)
        f = symfun(realSphericalHarmonics{k}, [mu, phi]);
        intRef = int(int(f, mu, -1, 1), phi, 0, 2 * pi);
        fAtPoints = double(f(muQuad, phiQuad));
        intQuad = sum(fAtPoints .* weights);
        assert(abs(intQuad - double(intRef)) < 1e-14);
    end
end

%% quadrature on upper half sphere
addpath('../basisFunctions/')
mu = sym('mu');
phi = sym('phi');
Nmax = 5;
for N = 2 : Nmax  % no rule for available for N = 1
    realSphericalHarmonics = getRealSphericalHarmonicsSym(N);
    [weights, muQuad, phiQuad] = sphericalQuadratureUpperHalf(N);
    for k = length(realSphericalHarmonics) - 2 * N : length(realSphericalHarmonics)
        f = symfun(realSphericalHarmonics{k}, [mu, phi]);
        intRef = int(int(f, mu, 0, 1), phi, 0, 2 * pi);
        fAtPoints = double(f(muQuad, phiQuad));
        intQuad = sum(fAtPoints .* weights);
        assert(abs(intQuad - double(intRef)) < 1e-14);
    end
end

%% quadrature on x > 0 sphere via rotation
addpath('../basisFunctions/')
outerNormalVector = [-1; 0; 0];
mu = sym('mu');
phi = sym('phi');
Nmax = 5;
for N = 2 : Nmax  % no rule for available for N = 1
    realSphericalHarmonics = getRealSphericalHarmonicsSym(N);
    [weights, muQuad, phiQuad] = sphericalQuadratureHalf(outerNormalVector, N);
    for k = length(realSphericalHarmonics) - 2 * N : length(realSphericalHarmonics)
        f = symfun(realSphericalHarmonics{k}, [mu, phi]);
        intRef = int(int(f, mu, -1, 1), phi, -pi / 2, pi/2);
        fAtPoints = double(f(muQuad, phiQuad));
        intQuad = sum(fAtPoints .* weights);
        assert(abs(intQuad - double(intRef)) < 1e-14);
    end
end

%% quadrature on y < 0 sphere via rotation
addpath('../basisFunctions/')
outerNormalVector = [0; 1; 0];
mu = sym('mu');
phi = sym('phi');
Nmax = 5;
for N = 2 : Nmax  % no rule available for N = 1
    realSphericalHarmonics = getRealSphericalHarmonicsSym(N);
    [weights, muQuad, phiQuad] = sphericalQuadratureHalf(outerNormalVector, N);
    for k = length(realSphericalHarmonics) - 2 * N : length(realSphericalHarmonics)
       f = symfun(realSphericalHarmonics{k}, [mu, phi]);
        intRef = int(int(f, mu, -1, 1), phi, pi, 2 * pi);
        fAtPoints = double(f(muQuad, phiQuad));
        intQuad = sum(fAtPoints .* weights);
        assert(abs(intQuad - double(intRef)) < 1e-14);
    end
end