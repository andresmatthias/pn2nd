function u = generateSymMomentVectorDerivative(spatialVar, N, spatialDimension)
% generate symbolic moment vector, indexed linearly starting from 0.
addpath('./basisFunctions/')
if spatialDimension == 1
    linIdxMax = N;
elseif (spatialDimension == 2) || (spatialDimension == 3)
    linIdxMax = DegOrder2linearIdx(N, N, spatialDimension);
else
    error('Invalid dimension!')
end

u = sym(zeros(linIdxMax + 1, 1));
for l = 0 : linIdxMax
    for s = spatialVar
        u(l + 1) = str2sym(sprintf('u%d_%s', l, s));
    end
end

end