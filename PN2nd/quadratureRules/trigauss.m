function tw=trigauss(n,alpha,beta)
% Downloaded from:
% https://www.math.unipd.it/~marcov/mysoft/subp/trigauss.m
% n: maximum trigonometric degree which is integrated exactly, i.e., 
% functions in f in span{1, cos(k*phi), sin(k*phi), 1<=k<=n} are integrated 
% exactly.
% Why called degree? c0 + c1*z + ... + cn*z^n is a complex polynomial of
% degree n; for z = e^{i*k*phi} with -n <= k <= n (restriction on unit 
% circle) we obtain the above representation in terms of cos(kx), sin(kx).
% 
%   For a demonstration on how to use this function,
%   see also SPHERICALQUADRATURETEST

% by Gaspare Da Fies, Alvise Sommariva and Marco Vianello
% University of Padova
% February 28, 2013

% computes the n+1 angles and weights of a trigonometric gaussian
% quadrature formula on [alpha,beta], 0<beta-alpha<=2*pi

% improves a previous version by G. Da Fies and M. Vianello

% input:
% n: trigonometric degree of exactness
% [alpha,beta]: angular interval, 0<beta-alpha<=2*pi

% output:
% tw: (n+1) x 2 array of (angles,weights)

% the formula integrates the canonical trigonometric basis with accuracy 
% from about 10^(-15) (small omega) to about 10^(-13) (omega-->pi) 
% up to n=300 


n=n+1;
omega=(beta-alpha)/2;
ab = r_subchebyshev(n,omega);
xw_symm_eigw = SymmMw(n,ab);
tw=trigauss_conversion(xw_symm_eigw,omega);
tw(:,1)=tw(:,1)+(beta+alpha)/2;




function ab=r_subchebyshev(n,omega)
%% Authors: 
%% Gerard Meurant and Alvise Sommariva 
%%
%% Date: June 3, 2012
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% INPUT:
% n    : NUMBER OF POINTS. 
% omega: ARC ANGLE.
%
% OUTPUT:
% ab   : COEFFS of THREE TERM RECURSION.
%
%--------------------------------------------------------------------------


N = n;
n = n - 1;


% modified Chebyshev moments by recurrence

if rem(N,2) == 1
    NN=N+1; nn=n+1;
else
    NN=N; nn=n;
end

mom=fast_moments_computation(omega,2*nn+1);

% recurrence coeffs of the monic Chebyshev polynomials
abm(:,1)=zeros(2*nn+1,1);
abm(:,2)=0.25*ones(2*nn+1,1); abm(1,2)=pi; abm(2,2)=0.5;

% recurrence coeffs for the monic OPS w.r.t. the weight function
% w(x)=2*sin(omega/2)/sqrt(1-sin^2(omega/2)*x^2)
% by the modified Chebyshev algorithm

% ab = chebyshev(NN+1,mom,abm);
ab = fast_chebyshev(NN,mom,abm);




function x = tridisolve(a,b,c,d)
%   TRIDISOLVE  Solve tridiagonal system of equations.
% From Cleve Moler's Matlab suite
% http://www.mathworks.it/moler/ncmfilelist.html

%     x = TRIDISOLVE(a,b,c,d) solves the system of linear equations
%     b(1)*x(1) + c(1)*x(2) = d(1),
%     a(j-1)*x(j-1) + b(j)*x(j) + c(j)*x(j+1) = d(j), j = 2:n-1,
%     a(n-1)*x(n-1) + b(n)*x(n) = d(n).
%
%   The algorithm does not use pivoting, so the results might
%   be inaccurate if abs(b) is much smaller than abs(a)+abs(c).
%   More robust, but slower, alternatives with pivoting are:
%     x = T\d where T = diag(a,-1) + diag(b,0) + diag(c,1)
%     x = S\d where S = spdiags([[a; 0] b [0; c]],[-1 0 1],n,n)

% optimized version 
x = d;
n = length(x);
bi = zeros(n,1);

for j = 1:n-1
 bi(j) = 1 / b(j);
 mu = a(j) * bi(j);
 b(j+1) = b(j+1) - mu * c(j);
 x(j+1) = x(j+1) - mu * x(j);
end

x(n) = x(n) / b(n);
for j = n-1:-1:1
 x(j) = (x(j) - c(j) * x(j+1)) * bi(j);
end




function ab=fast_chebyshev(N,mom,abm)
%SUBP_MOD_CHEBYSHEV Modified Chebyshev algorithm
% this works only for the subperiodic weight function
%
% From Gautschi's code (simplified)
% Mar 2012
%

ab = zeros(N,2);
sig = zeros(N+1,2*N);

ab(1,2) = mom(1);

sig(1,1:2*N) = 0; 
sig(2,:) = mom(1:2*N);

for n = 3:N+1
  for m = n-1:2*N-n+2
     sig(n,m) = sig(n-1,m+1) + abm(m,2) * sig(n-1,m-1) - ab(n-2,2) * sig(n-2,m);
  end
 
  ab(n-1,2) = sig(n,n-1) / sig(n-1,n-2);
end




function mom=fast_moments_computation(omega,n)

mom=zeros(1,n+1);
mom(1)=2*omega; % FIRST MOMENT.

if(n>=2)

    if(omega<=1/4*pi)
        l=10;
    elseif(omega<=1/2*pi)
        l=20;
    elseif(omega<=3/4*pi)
        l=40;
    else
        if omega == pi
            l=2*ceil(10*pi);
        else
            l=2*ceil(10*pi/(pi-omega));
        end
    end

    
    temp=(2:2:n+2*l-2); % AUXILIAR VECTORS.
    temp2=temp.^2-1;

    dl=1/4 -1./(4*(temp-1)); % DIAGONALS.
    dc=1/2 -1/sin(omega/2)^2 -1./(2*temp2);
    du=1/4 +1./(4*(temp+1));

    d=4*cos(omega/2)/sin(omega/2)./temp2'; % COMPUTING KNOWN TERM.
    d(end)=d(end);                         % PUT LAST MOMENT NULL.

    z=tridisolve(dl(2:end),dc,du(1:end-1),d); % SOLVE SYSTEM.
    mom(3:2:n+1)=z(1:floor(n/2)); % SET ODD MOMENTS.

end

mom=mom';

normalized = 0;

if normalized == 0
    M=length(mom);
    kk=2.^(-((1:2:M)-2))'; kk(1)=1;
    v=ones(M,1);
    v(1:2:M)=kk;
    mom=v.*mom;
end




function xw=SymmMw(N,ab)
%SymmMw computation of the nodes and weights for a symmetric weight
%function
% this version uses the reduced matrix and eig and
% computation of weights with the 3-term recurrence

%
% see: Fast variants of the Golub and Welsch algorithm for symmetric
% weight functions by G. Meurant and A. Sommariva (2012)

% Input
% N : cardinality of the rule
% ab: 3-term recurrence for the orthogonal polynomials
% same as in OPQ
% ab(1,2) is the 0th moment

% Output
% xw : xw(:,1) nodes, xw(:,2) weights of the quadrature rule
% nloop: number of iterations in QR

%
% Authors G. Meurant and A. Sommariva
% June 2012
%

N0 = size(ab,1);
if N0 < N
 error('SymmMw: input array ab is too short')
end

na = norm(ab(:,1));
if na > 0
 error('SymmMw: the weight function must be symmetric')
end

% computation of the reduced matrix in vectors (a,b)

if mod(N,2) == 0
    even = 1;
    Nc = N / 2;
else
    even = 0;
    Nc = fix(N / 2) +1;
end


absd = ab(:,2);
absq = sqrt(absd);

a = zeros(1,Nc);
b = a;

switch even
    case 1
        % N even
        a(1) = absd(2);
        b(1) = absq(2) * absq(3);
 
        k = (2:Nc-1);
        a(k) = absd(2*k-1) + absd(2*k);
        b(k) = absq(2*k) .* absq(2*k+1);
        a(Nc) = absd(N) + absd(N-1);
        start = 1;
        
       J = diag(a) + diag(b(1:Nc-1),1) + diag(b(1:Nc-1),-1);
       t = sort(eig(J));
       w = weights_3t(t',a,b);
       % w are the squares of the first components
       w = w' / 2;
    case 0
        % N odd
        a(1) = absd(2);
        b(1) = absq(2) * absq(3);
        
        k = (2:Nc-1);
        a(k) = absd(2*k-1) + absd(2*k);
        b(k) = absq(2*k) .* absq(2*k+1);
        a(Nc) = absd(N);
        start = 2;

        % the first node must be zero
        J = diag(a) + diag(b(1:Nc-1),1) + diag(b(1:Nc-1),-1);
        t = sort(eig(J));
        t(1) = 0;
        w = weights_3t(t',a,b);
        w = [w(1); w(2:end)' / 2];
    otherwise
        error('this is not possible')
end

xwp = sqrt(t);

xw(:,1) = [-xwp(end:-1:start,1); xwp];
xw(:,2) = ab(1,2) * ([w(end:-1:start); w]);




function tw=trigauss_conversion(xw,omega)

tw(:,1)=2*asin(sin(omega/2)*xw(:,1));
tw(:,2)=xw(:,2);




function w=weights_3t(t,a,b)
%WEIGHTS_3T squares of the 1st components of eigenvectors from the 3-term  
% recurrence relation of the orthogonal polynomials
%

% Input
% t: nodes
% a,b coefficients of the 3-term recurrence
%
% Ouput
% w: squares of the first components of the eigenvectors
%

%
% Authors G. Meurant and A. Sommariva
% June 2012
%

N = length(t);

P = zeros(N,N);
P(1,:) = ones(1,N);
P(2,:) = (t - a(1)) / b(1);

for k = 3:N
 k1 = k - 1;
 k2 = k - 2;
 P(k,:) = ((t - a(k1)) .* P(k1,:) - b(k2) * P(k2,:)) / b(k1);
end

P2 = P .* P;

w = 1 ./ sum(P2);


