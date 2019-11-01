addpath(genpath('../PN2nd'));
clear
close all
[weights, mu, phi] = sphericalQuadratureFull(20);
[realSphericalHarmonicsAtPoints] = evalRealSphericalHarmonics([mu;phi], 2);
Omega = sphereParam2Cartesian(mu, phi);
K = convhull(Omega');

k = 2;
X = reshape(Omega(1, K(:)), size(K, 1), 3);
Y = reshape(Omega(2, K(:)), size(K, 1), 3);
Z = reshape(Omega(3, K(:)), size(K, 1), 3);
f = reshape(realSphericalHarmonicsAtPoints(k, K(:)), size(K, 1), 3);
f = mean(f, 2);


patch(X', Y', Z', f')
axis equal vis3d
view(3)