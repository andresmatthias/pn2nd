function R = getGeomRotationa2b(a, b)
% GETGEOMROTATIONA2B Get rotation matrix, that transforms a given unit  
%   vector a to a desired unit vector b by geometric rotation, i.e.,
%   b = R * a with det(R) = 1, R' * R = I.
% 
%   For a demonstration on how to use this function,
%   see also MISCTEST
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

if (abs(norm(a) - 1) + abs(norm(b) - 1)) > 1e-14
   error('Vectors need to be normalized!') 
end

if norm(a + b) < 1e-15
    %formula is not applicable in this case but integrals can be
    %computed analogously
    R = diag([1, 1, -1]);
else
    v = cross(a, b);
    c = dot(a, b);
    vx = [0, -v(3), v(2);
          v(3),0, -v(1);
          -v(2), v(1), 0];
    R = eye(3) + vx + vx * vx / (1 + c);
end    
end