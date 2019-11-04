% DRIFTKERNEL Check the projected scattering operator for the special
%   kernel, that would yield a drift term in the second-order formulation.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

clear
addpath('../PN2nd/basisFunctions/')
addpath('../tools/')
syms mu phi mup phip

k =    sym(3) / sym(1900) / sym(pi) * (- 15 * mu - 15 * mup ...
                        - 75 * mu^2 * mup^2 - 27 * mu * mup ...
                        + 45 * mu * mup ^2 + 45 * mu^2 * mup ...
                        + 25 * mu^2 + 25 * mup^2 + sym(150));

%%
I = int(int(k, mu, -1, 1), phi, 0, 2*pi);
fprintf('Should be 1:\n')
disp(I)

%%

[MU, MUP] = meshgrid(linspace(-1, 1, 100), linspace(-1, 1, 100));
kFun = matlabFunction(k);
Z = kFun(MU, MUP);

cmap = viridis(256);
colormap(cmap);
surf(MU, MUP, Z); 
view(2); shading interp
colorbar
                         


%%
b = getRealSphericalHarmonicsSym(3);
idxEven = linearIdxOfEvenBasis(3, 3) + 1;
idxOdd = linearIdxOfOddBasis(3, 3) + 1;
be = b(idxEven);
bo = b(idxOdd);

for j = 1 : length(idxOdd)
    bo{j} = subs(bo{j}, [mu, phi], [mup, phip]);
end

Ceo = sym(zeros(length(idxEven), length(idxOdd)));
nzFlag = 0;
for i = 1 : length(idxEven)
    for j = 1 : length(idxOdd)
        fprintf('i/j: %d/%d\n', i, j)
        tmp = int(int(be{i} * bo{j} * k, mu, -1, 1), phi, 0, 2*pi);
        Ceo(i, j) = simplify(int(int(tmp, mup, -1, 1), phip, 0, 2*pi));
        disp(Ceo(i, j))
        if vpa(Ceo(i,j)) ~=0
            nzFlag = 1;
        end
    end
end

if nzFlag ~=0
   fprintf('non-zero entry found\n'); 
else
   fprintf('all entries zero\n')
end
disp(Ceo)




