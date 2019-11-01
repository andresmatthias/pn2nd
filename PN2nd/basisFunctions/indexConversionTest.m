%INDEXCONVERSIONTEST Perform unit tests regarding the conversion of linear
%   basis indices (when serially numbering all basis functions, depending 
%   on symmetry assumptions (spatialDimension)) into corresponding
%   degree-order tuples and vice versa.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

%% conversion from DegreeOrder to linear index and back
lmax = 50;
linIdxTrue = 0 : lmax;
for spatialDimension = [1, 2, 3]
    [l, m] = linearIdx2DegOrder(linIdxTrue, spatialDimension);
    linIdx = degOrder2linearIdx(l, m, spatialDimension);
    assert(isequal(linIdx, linIdxTrue))
end

%% conversion from DegreeOrder index in 1D to linear index 
% (only those spherical harmonics with m == 0
% contained in basis set)
spatialDimension = 1;
lmax = 50;
lList = 0 : lmax;
mList = zeros(size(lList));
linIdxConv = degOrder2linearIdx(lList, mList, spatialDimension);
linIdxTrue = lList;
assert(isequal(linIdxConv, linIdxTrue));

%% conversion from linear idx to DegreeOrder index in 1D
spatialDimension = 1;
linIdx = 0 : 50;
[lConv, mConv] = linearIdx2DegOrder(linIdx, spatialDimension);
lListTrue = linIdx;
mListTrue = zeros(size(lListTrue));
assert(isequal(lListTrue, lConv));
assert(isequal(mListTrue, mConv));

%% conversion from DegreeOrder index in 2D to linear index 
% (only those spherical harmonics with mod(l+m, 2) == 0 (even parity)
% contained in basis set)
spatialDimension = 2;
lmax = 50;
linIdx = [];
lList = [];
mList = [];
cnt = 0;
for l = 0:lmax
    for m = -l : 2 : l % skip basisfunctions with odd parity
        lList(end + 1) = l;
        mList(end + 1) = m;
        linIdx(end + 1) = cnt;
        cnt = cnt + 1;
    end
end
linIdxConv = degOrder2linearIdx(lList, mList, spatialDimension);
assert(isequal(linIdxConv, linIdx));

%% conversion from linear idx to DegreeOrder index in 2D
spatialDimension = 2;
lmax = 50;
linIdx = [];
lList = [];
mList = [];
cnt = 0;
for l = 0:lmax
    for m = -l : 2 : l % skip basisfunctions with odd parity
        lList(end + 1) = l;
        mList(end + 1) = m;
        linIdx(end + 1) = cnt;
        cnt = cnt + 1;
    end
end
[lback, mback] = linearIdx2DegOrder(linIdx, spatialDimension);
assert(isequal(lback, lList));
assert(isequal(mback, mList));


%% conversion from DegreeOrder index in 3D to linear index
% all spherical harmonics contained in basis set
spatialDimension = 3;
lmax = 50;
linIdx = [];
lList = [];
mList = [];
cnt = 0;
for l = 0:lmax
    for m = -l : 1 : l
        lList(end + 1) = l;
        mList(end + 1) = m;
        linIdx(end + 1) = cnt;
        cnt = cnt + 1;
    end
end
linIdxConv = degOrder2linearIdx(lList, mList, spatialDimension);
assert(isequal(linIdxConv, linIdx));

%% conversion from linear idx to DegreeOrder index in 3D
spatialDimension = 3;
lmax = 50;
linIdx = [];
lList = [];
mList = [];
cnt = 0;
for l = 0:lmax
    for m = -l : 1 : l % skip basisfunctions with odd parity
        lList(end + 1) = l;
        mList(end + 1) = m;
        linIdx(end + 1) = cnt;
        cnt = cnt + 1;
    end
end
[lback, mback] = linearIdx2DegOrder(linIdx, spatialDimension);
assert(isequal(lback, lList));
assert(isequal(mback, mList));
