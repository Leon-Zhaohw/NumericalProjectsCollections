function [ distanceField ] = fluidNodesToDistanceField(fluidNodes)
%FLUIDNODESTODISTANCEFIELD Convert a boolean matrix denoting fluid nodes 
%       to a distance field.
%   return distanceField: double matrix(yNodes, xNodes)
%   fluidNodes: boolean matrix(yNodes, xNodes)

    % calculate distanceField from waterbody
    distanceField = zeros(size(fluidNodes));
    distancesToAir = (bwdist(~fluidNodes) - 0.5);
    distancesToWater = (bwdist(fluidNodes) - 0.5) .* -1;
    distanceField(~fluidNodes) = distancesToWater(~fluidNodes);
    distanceField(fluidNodes) = distancesToAir(fluidNodes);

end

