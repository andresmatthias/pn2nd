function [outgoingOrdinateIdx] = getOutgoingOrdinateIdx1D(...
    ingoingOrdinateIdx, ordinates)
% GETOUTGOINGORDINATEIDX1D Get the index of the reflected discrete
%   ordinate in 1D. Assuming that ordinates are symmetrical around 
%   mu = 0, we can look for the closest correspondence and don't need to 
%   interpolate.
%
% For details, see our publication on arXiv:
% The second-order formulation of the PN equations with Marshak boundary conditions
% by Matthias Andres and Florian Schneider
% 1 Nov 2019
% https://arxiv.org/abs/1911.00468
%

ingoingOrdinatesMu = ordinates.mu(ingoingOrdinateIdx);
ingoingOrdinatesPhi = ordinates.phi(ingoingOrdinateIdx);
if any(ingoingOrdinateIdx) % in PN_FV this array might be empty
    outgoingOrdinatesMu = -ingoingOrdinatesMu; % 1D !!!

     %find closest correspondence in discrete_ordinates     
     aux = outgoingOrdinatesMu(:) - ordinates.mu(:)';
     auxPhi = ingoingOrdinatesPhi(:) - ordinates.phi(:)';
     
     [m, outgoingOrdinateIdx] = min(abs(aux) + abs(auxPhi), [], 2);
     outgoingOrdinateIdx = outgoingOrdinateIdx';
     if any(m>1e-15)
         error('sth went wrong with reflection!')
     end
else
    outgoingOrdinateIdx  = zeros(1, 0);
end
end

