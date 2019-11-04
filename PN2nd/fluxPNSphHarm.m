function [Tx, Ty, Tz] = fluxPNSphHarm(N, spatialDimension)
% FLUXPNSPHHARM Compute flux matrices of the PN system with the real
%   spherical harmonics as basis functions.
%
%   (Flux matrices: Mx * dxu + My * dyu + Mz * dzu = ...)
%   Assemble full system and reduce dimension according to symmetry
%   assumptions (spatial dimension) at the end.
%
%   Implementation based on:
%       StaRMAP---A Second Order Staggered Grid Method for Spherical
%       Harmonics Moment Equations of Radiative Transfer; Benjamin
%       Seibold, Martin Frank; 2014; ACM Transactions on Mathematical
%       Software; Volume 41; 
%       doi: 10.1145/2590808
%       
%   For a demonstration on how to use this function,
%   see also FLUXMATRIXTEST
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

% addpath('./basisFunctions/')
spatialDimensionFull = 3;
nMomentsFull = getNumberOfBasisFunctions(N, spatialDimensionFull);

[lListFull, mListFull] = linearIdx2DegOrder(0 : nMomentsFull - 1, spatialDimensionFull); 
degOrder = [lListFull', mListFull'];

Mxc = computeMxComplexFull();
Myc = computeMyComplexFull();
Mzc = computeMzComplexFull(); 
%% transform the system matrices to real system
S = trafoComp2RealBasis(nMomentsFull);
Tx = conj(S) * Mxc * transpose(S);
Ty = conj(S) * Myc * transpose(S);
Tz = conj(S) * Mzc * transpose(S);

if (norm(imag(Tx), 'inf') + norm(imag(Ty), 'inf') + norm(imag(Tz), 'inf')) > 1e-15
   error('Something wrong with flux matrices: they should be real!')
else
    Tx = real(Tx);
    Ty = real(Ty); 
    Tz = real(Tz);
end

%% reduction of dimension due to symmetry
idxReduced = getReducedIdx(lListFull, mListFull, spatialDimension); 
Tx = Tx(idxReduced, idxReduced);
Ty = Ty(idxReduced, idxReduced);
Tz = Tz(idxReduced, idxReduced);


%%
function Mxc = computeMxComplexFull()
    Mxc = zeros(nMomentsFull, nMomentsFull);
    %% equation (5) in StaRMAP (see above)
    for i=1 : nMomentsFull
        l = degOrder(i, 1); m = degOrder(i, 2);
        if (l > 0) && (abs(m - 1) <= (l - 1))
            j1 = degOrder2linearIdx(l - 1, m - 1, spatialDimensionFull) + 1; % Matlab indexing + 1
            Mxc(i, j1) = -c(l - 1, m - 1);
        end
        if (l < N) && (abs(m - 1) <= l + 1)
            j2 = degOrder2linearIdx(l + 1, m - 1, spatialDimensionFull) + 1;
            Mxc(i, j2) = d(l + 1, m - 1);
        end
        if (l > 0) && (abs(m + 1) <= l - 1)
            j3 = degOrder2linearIdx(l - 1, m + 1, spatialDimensionFull) + 1;
            Mxc(i, j3) = e(l - 1, m + 1);
        end
        if (l < N) && (abs(m + 1) <= l + 1)
           j4 = degOrder2linearIdx(l + 1, m + 1, spatialDimensionFull) + 1;
            Mxc(i, j4) = -f(l + 1, m + 1);
        end
    end
    Mxc = 1 / 2 * Mxc;
end

function Myc = computeMyComplexFull()
    Myc =zeros(nMomentsFull, nMomentsFull);
    %% equation (5) in StaRMAP (see above)
    for i=1 : nMomentsFull
        l = degOrder(i, 1); m = degOrder(i, 2);
        if (l > 0) && (abs(m - 1) <= (l - 1))
            j1 = degOrder2linearIdx(l - 1, m - 1, spatialDimensionFull) + 1; % Matlab indexing + 1
            Myc(i, j1) = c(l - 1,m - 1);
        end
        if (l < N) && (abs(m - 1) <= l + 1)
            j2 = degOrder2linearIdx(l + 1, m - 1, spatialDimensionFull) + 1;
            Myc(i, j2) = -d(l + 1, m - 1);
        end
        if (l > 0) && (abs(m + 1) <= l - 1)
            j3 = degOrder2linearIdx(l - 1, m + 1, spatialDimensionFull) + 1;
            Myc(i, j3) = e(l - 1, m + 1);
        end
        if (l < N) && (abs(m + 1) <= l + 1)
           j4 = degOrder2linearIdx(l + 1, m + 1, spatialDimensionFull) + 1;
            Myc(i, j4) = -f(l + 1, m + 1); 
        end
    end
    Myc = 1i / 2 * Myc;
end

function Mzc = computeMzComplexFull()
    Mzc =zeros(nMomentsFull, nMomentsFull);
    %% equation (5) in https://dl.acm.org/citation.cfm?id=2590808
    for i=1 : nMomentsFull
        l = degOrder(i, 1); m = degOrder(i, 2);
        if (l > 0) && (abs(m) <= l - 1)
            j5 = degOrder2linearIdx(l - 1, m, spatialDimensionFull) + 1;
            Mzc(i, j5) = a(l - 1, m);
        end
        if l < N  && (abs(m) <= l + 1)
            j6 = degOrder2linearIdx(l + 1, m, spatialDimensionFull) + 1;
            Mzc(i, j6) = b(l + 1, m);
        end
    end
end

function S = trafoComp2RealBasis(nMoments)
%% compute transformation matrix S to convert complex into real spherical harmonics
% S_l^m = (-1)^m/sqrt(2) * (Y_l^m + (-1)^m Y_l^{-m})               if m > 0
%         Y_l^0                                                    if m = 0
%         -i * (-1)^m/sqrt(2) * (Y_l^{|m|} - (-1)^m Y_l^{-|m|})    if m < 0
    S = zeros(nMoments, nMoments);
    for i = 1 : nMoments
        l = degOrder(i, 1); m = degOrder(i, 2);
        if m == 0
            S(i, i) = 1; 
        elseif m > 0
            S(i, i) = (-1)^m / sqrt(2); % prefactor of Y_l^m
            idx = degOrder2linearIdx(l, -m, spatialDimensionFull) + 1;
            S(i, idx) = 1 /  sqrt(2); % prefactor of Y_l^{-m}
        elseif m < 0
            S(i, i) = 1i / sqrt(2); % prefactor of Y_l^{-|m|}
            idx = degOrder2linearIdx(l, abs(m), spatialDimensionFull) + 1;
            S(i, idx) = -1i * (-1)^abs(m) / sqrt(2); % prefactor of Y_l^{|m|}
        end
    end
end

end
function y = a(l, m)
    y = sqrt( (l - m + 1) * (l + m + 1) / (2 * l + 3) / (2 * l + 1));
end

function y = b(l, m)
    y = sqrt( (l - m) * (l + m) / (2 * l + 1) / (2 * l - 1) );
end

function y = c(l, m)
    y = sqrt( (l + m + 1) * (l + m + 2) / (2 * l + 3) / (2 * l + 1) );
end

function y = d(l, m)
    y = sqrt( (l - m) * (l - m - 1) / (2 * l + 1) / (2 * l - 1) );
end

function y = e(l, m)
    y = sqrt( (l - m + 1) * (l - m + 2) / (2 * l + 3) / (2 * l + 1) );
end

function y = f(l, m)
    y = sqrt( (l + m) * (l + m - 1) / (2 * l + 1) / (2 * l - 1));
end
