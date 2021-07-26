function [ correctedVelocity ] = setAirVelocity( ...
    velocity, distanceField, obstacle )
%SETAIRVELOCITY This function correct velocities in the air. It copies
%       the value from the nearest fluid cell. The path to the fluid cell 
%       cannot go through a solid node.
%   return correctedVelocity: matrix(yNodes, xNodes, velocityInDim)
%   velocity: matrix(yNodes, xNodes, velocityInDim)
%   distanceField: matrix(yNodes, xNodes)
%   obstacle: boolean matrix(yNodes, xNodes)

    ySize = size(obstacle, 1);
    xSize = size(obstacle, 2);

    waterBody = distanceFieldToFluidNodes(distanceField);
    
    % infinity here is a distance bigger than the maximum possible 
    % distance through normal nodes
    infiniteDist = xSize * ySize + 1000;
    
    weights = ones(size(distanceField));
    weights(obstacle) = infiniteDist;
    
    % get nearest water cell indieces
    % use graydist instead of bwdist for weighted nodes (obstacles)
    [distToFluid, nearestFluidCellIndieces] = bwdist(waterBody);
    weightedDistToFluid = graydist(weights, waterBody);
    
    correctedVelocity = velocity;
    
    % no air left, return input velocity
    if max(max(max(nearestFluidCellIndieces))) == 0
        return;
    end
    
    velX = velocity(:,:,1);
    vely = velocity(:,:,2);

    % get values of nearest water cells
    correctedAirVelocity = zeros(size(velocity));
    correctedAirVelocity(:,:,1) = velX(nearestFluidCellIndieces);
    correctedAirVelocity(:,:,2) = vely(nearestFluidCellIndieces);
    
    % set velocity in air nodes behind obstacles to zero
    blockedAirNodes = repmat(weightedDistToFluid >= infiniteDist, [1 1 2]);
    correctedAirVelocity(blockedAirNodes) = 0;

    % set velocity in air nodes in obstacles to zero
    correctedAirVelocity(repmat(obstacle, [1 1 2])) = 0;
    
    % copy values of nearest water cells to air cells
    inflatedWaterBody = repmat(waterBody, [1 1 2]);
    correctedVelocity(~inflatedWaterBody) = ...
    correctedAirVelocity(~inflatedWaterBody);

end

