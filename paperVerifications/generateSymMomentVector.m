function u = generateSymMomentVector(N, dimension)
% generate symbolic moment vector, indexed linearly starting from 0.
if dimension == 1
    linIdxMax = N;
elseif (dimension == 2) || (dimension == 3)
    linIdxMax = DegOrder2linearIdx(N, N, dimension);
else
    error('Invalid dimension!')
end

u = sym(zeros(linIdxMax + 1, 1));
for l = 0 : linIdxMax
    u(l + 1) = str2sym(sprintf('u%d', l));
end

end