function [ newVelocityField ] = enforceBoundaryCondition( velocityField, obstacle )
%ENFORCEBOUNDARYCONDITION This function ensures that the velocity in 
%       normal direction from the solids are zero. This is the free-slip
%       boundary condition.
%   return newVelocityField: matrix(yNodes, xNodes, dim)
%   velocityField: matrix(yNodes, xNodes, dim)
%   obstacle: boolean matrix(yNodes, xNodes)

    mSize = size(velocityField, 1);
    nSize = size(velocityField, 2);
    
    newVelocityField = velocityField;
    
    paddedObstacle = padarray(obstacle, [1 1]);
    
    % Find neighboring obstacles
    obstacleLeft = paddedObstacle(2 : mSize + 1, 3 : nSize + 2);
    obstacleRight = paddedObstacle(2 : mSize + 1, 1 : nSize);
    obstacleUp = paddedObstacle(3 : mSize + 2, 2 : nSize + 1);
    obstacleDown = paddedObstacle(1 : mSize, 2 : nSize + 1);
    
    velocityToDelete = false(size(velocityField));
    
    % free-slip condition: velocity in normal direction must be zero
    velocityToDelete(:,:,1) = obstacleLeft | obstacleRight;
    velocityToDelete(:,:,2) = obstacleDown | obstacleUp;
    
    newVelocityField(velocityToDelete) = 0;
    
end

