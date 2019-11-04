function edgesPerElement = getEdgesPerElement(connectivityList, edgeList)
% GETEDGESPERELEMENT Get the edge IDs per element.
%  
%   For each element, the k-th (k=1,2,3) edge is opposite of the k-th point
%   in the connectivity list. 
% 
%   For a demonstration on how to use this function,
%   see also MESHAUXTEST
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%


%% adjacent edges
% ordered: opposite of nodes in connectivity_list
% [A,B,C] = connectivityList(k, :);
% element_edges : [BC, CA, AB]
fprintf('\nGet edges ...\n')
[~, edgesAB] = ismember(sort(connectivityList([1,2], :), 1)', edgeList', 'rows');
[~, edgesBC] = ismember(sort(connectivityList([2,3], :), 1)', edgeList', 'rows');
[~, edgesCA] = ismember(sort(connectivityList([3,1], :), 1)', edgeList', 'rows');
edgesPerElement = [edgesBC'; edgesCA'; edgesAB'];
end