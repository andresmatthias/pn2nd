% MISCDISCRETEORDINATESTEST Perform unit tests regarding the discrete
% ordiantes method.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

%% index conversion: linear to tuple
nOrdinates = 7;
linIdx = randi(1000, 1, 100);
tupIdx = lin2tupDO(linIdx, nOrdinates);
linIdxRet = tup2linDO(tupIdx, nOrdinates);
assert(all(linIdxRet == linIdx))

%% index conversion: linear to tuple
nOrdinates = 7;
tupIdx = [randi(1000, 1, 100);
          randi(nOrdinates, 1, 100)];
linIdx = tup2linDO(tupIdx, nOrdinates);
tupIdxRet = lin2tupDO(linIdx, nOrdinates);
assert(all(tupIdxRet(:) == tupIdx(:)))

%% weights sum up to 4 pi
addpath('../PN2nd/quadratureRules/')
addpath('../PN2nd/')
maxExactDegree = 10;
kernelName = 'isotropic';
for dimension = [1, 2, 3]
    ordinates = getOrdinates(maxExactDegree, dimension, kernelName);
    assert(abs(sum(ordinates.weights) - 4 * pi) < 1e-13);
end

