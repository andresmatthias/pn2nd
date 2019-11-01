function [outgoingOrdinateIdx, weights] = getOutgoingOrdinateIdxInterp2D(...
    ingoingOrdinateIdx, outerNormal, ordinates)
% GETOUTGOINGORDINATEIDXINTERP2D For a given outer normal vector at the
%   boundary and an ingoing direction (discrete ordinate), compute the 
%   original outgoing direction which was reflected at the boundary. 
%   Interpolate this outgoing direction based on a given set of discrete  
%   ordinates.
%
%See The second-order formulation of the P_N equations with Marshak boundary conditions 
%Matthias Andres, Florian Schneider
%

ingoingOrdinates_xy = ordinates.Cartesian(1:2, ingoingOrdinateIdx);
ingoingOrdinates_z = ordinates.Cartesian(3, ingoingOrdinateIdx);
if any(ingoingOrdinateIdx) 
    if nargin <= 3
        if any(dot(ingoingOrdinates_xy, outerNormal, 1) >= 0)
            error('Ordinate is not pointing inward!')
        end
    end
    outgoingOrdinates_xy = ingoingOrdinates_xy -...
        2 * dot(ingoingOrdinates_xy, outerNormal, 1) .* outerNormal;

    [outgoingOrdinateIdx, weights] = distanceBasedInterpolationUnitSphere(...
        ordinates.Cartesian, [outgoingOrdinates_xy; ingoingOrdinates_z]);
else
    outgoingOrdinateIdx  = zeros(3, 0);
    weights = nan(3, 0);
end
end

