function [kernel] = henyeyGreenstein(g)
% HENYEYGREENSTEIN Define the Henyey-Greenstein kernel function with
%   anisotropy factor g.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

% theta: angle between two ordinates
cos_theta = @(vx, vy, vz, wx, wy, wz) (vx .* wx + vy .* wy + vz .* wz) ./...
    (sqrt(vx.^2 + vy.^2 + vz.^2) .* sqrt(wx.^2 + wy.^2 + wz.^2));
kernel.fun = @(vx, vy, vz, wx, wy, wz) 1 / 4 / pi * (1 - g^2) ./ ...
    ((1 + g^2 - 2 * g * cos_theta(vx, vy, vz, wx, wy, wz)).^(3 / 2));
kernel.name = sprintf('henyeyGreenstein_g%1.3f', g);
kernel.name = strrep(kernel.name, '.', ',');
end

