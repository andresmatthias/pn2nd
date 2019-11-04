function [x] = choleskySolve(A, b)
% CHOLESKYSOLVE Solve A*x = b for symmetric A, with A = L*L.' and L lower 
%   triangular matrix.
% 
%   For a demonstration on how to use this function,
%   see also MISCTEST
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

n = size(A, 1);
if isnumeric(A)
    L = chol(A, 'lower');
else
    L = chol(A, 'lower', 'real', 'nocheck');
end

% Ly = r (lower triangle); L'x = y (upper triangle)
y = 0 * b + 0 * A(1,1); % in case A or b is symbolic
y(1, :) = b(1, :) ./ L(1,1);
for i = 2:n
   fprintf('forward: %s, %d / %d\n', datetime('now'), i, n)
   S = 0 * b(1,:);
   for j = 1:i-1
      S = S + L(i, j) .* y(j, :);
   end
   y(i, :) = (b(i, :) - S) ./ L(i, i);
   y(i, :) = y(i, :);
end

%
x = 0 * b + 0 * A(1,1); % in case A or b is symbolic
R = L.';
x(n, :) = y(n, :) ./ R(n, n);
for i = n-1: -1: 1
   fprintf('backward: %s, %d / %d\n', datetime('now'), i, n)
   S = 0 * y(1, :);
   for j = i+1:n
      S = S + R(i, j) .* x(j, :); 
   end
   x(i, :) = (y(i, :) - S) ./ R(i, i);
   x(i, :) = x(i, :);
end
end
