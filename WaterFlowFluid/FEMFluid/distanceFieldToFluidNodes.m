function [ fluidNodes ] = distanceFieldToFluidNodes(distanceField)
%DISTANCEFIELDTOFLUIDNODES Get nodes inside the water
%   return fluidNodes: boolean matrix(yNodes, xNodes)
%   distanceField: matrix(yNodes, xNodes)

    % get water filled areas
    fluidNodes = (distanceField >= 0);

end

