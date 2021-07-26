function [ divergence ] = calcDivergence( velocity, obstacle )
%CALCDIVERGENCE Calculate the divergence from a velocity field.
%   return divergence: matrix(yNodes, xNodes)
%   velocity: matrix(yNodes, xNodes, dim)
%   obstacle: boolean matrix(yNodes, xNodes)

    ySize = size(velocity, 1);
    xSize = size(velocity, 2);
    
    % make look-up at borders easier
    % same effect as clamped texture coordinates in opengl
    paddedVelocity = padarray(velocity, [1 1 0], 'replicate');
    paddedObstacle = padarray(obstacle, [1 1], 'replicate');
    
    % Find neighboring velocity
    velocityUp = paddedVelocity(2 : ySize + 1, 3 : xSize + 2, :);
    velocityDown = paddedVelocity(2 : ySize + 1, 1 : xSize, :);
    velocityRight = paddedVelocity(3 : ySize + 2, 2 : xSize + 1, :);
    velocityLeft = paddedVelocity(1 : ySize, 2 : xSize + 1, :);
    
    % Find neighboring obstacles
    obstacleUp = logical(paddedObstacle(2 : ySize + 1, 3 : xSize + 2));
    obstacleDown = logical(paddedObstacle(2 : ySize + 1, 1 : xSize));
    obstacleRight = logical(paddedObstacle(3 : ySize + 2, 2 : xSize + 1));
    obstacleLeft = logical(paddedObstacle(1 : ySize, 2 : xSize + 1));
    
    % Set velocities to 0 for solid cells
    velocityUp(obstacleUp) = 0;
    velocityDown(obstacleDown) = 0;
    velocityRight(obstacleRight) = 0;
    velocityLeft(obstacleLeft) = 0;
    
    % difference of x + difference of y around the cell
    divergence = 0.5 * (velocityRight(:,:,2) - velocityLeft(:,:,2) + ...
        velocityUp(:,:,1) - velocityDown(:,:,1));
    
end

