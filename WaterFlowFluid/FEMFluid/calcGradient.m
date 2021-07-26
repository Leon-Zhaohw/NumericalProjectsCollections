function [ gradient ] = calcGradient( pressure, obstacle )
%CALCGRADIENT Calculate the pressure gradient
% based on: https://github.com/candycat1992/2DFluidSim
%   return newVelocityField: matrix(yNodes, xNodes, dim)
%   pressure: matrix(yNodes, xNodes)
%   obstacle: boolean matrix(yNodes, xNodes)

    mSize = size(pressure, 1);
    nSize = size(pressure, 2);
    
    % make look-up at borders easier
    % same effect as clamped texture coordinates in opengl
    paddedPressure = padarray(pressure, [1 1]);
    paddedObstacle = padarray(obstacle, [1 1]);
    
    % Find neighboring pressure
    pressureLeft = paddedPressure(2 : mSize + 1, 3 : nSize + 2);
    pressureRight = paddedPressure(2 : mSize + 1, 1 : nSize);
    pressureUp = paddedPressure(3 : mSize + 2, 2 : nSize + 1);
    pressureDown = paddedPressure(1 : mSize, 2 : nSize + 1);
    
    % Find neighboring obstacles
    obstacleLeft = logical(paddedObstacle(2 : mSize + 1, 3 : nSize + 2));
    obstacleRight = logical(paddedObstacle(2 : mSize + 1, 1 : nSize));
    obstacleUp = logical(paddedObstacle(3 : mSize + 2, 2 : nSize + 1));
    obstacleDown = logical(paddedObstacle(1 : mSize, 2 : nSize + 1));
    
    % Use center pressure for solid cells
    pressureLeft(obstacleLeft) = pressure(obstacleLeft);
    pressureRight(obstacleRight) = pressure(obstacleRight);
    pressureUp(obstacleUp) = pressure(obstacleUp);
    pressureDown(obstacleDown) = pressure(obstacleDown);
    
    % Note where forward / backward difference is used
    % divide by 2 only for central differences
    differenceFactorHorizontal = 2 .* ones(size(pressure));
    differenceFactorHorizontal(obstacleLeft | obstacleRight) = 1;
    differenceFactorVertical = 2 .* ones(size(pressure));
    differenceFactorVertical(obstacleUp | obstacleDown) = 1;

    % calc gradient with central differences
    gradient = zeros(mSize, nSize, 2);
    gradient(:, :, 1) = (pressureRight(:, :) - pressureLeft(:, :)) ./ ...
        differenceFactorHorizontal; 
    gradient(:, :, 2) = (pressureDown(:, :) - pressureUp(:, :)) ./ ...
        differenceFactorVertical;
    
    % set gradient in obstacles to 0
    % velocity in obstacles will be set to 0 anyway, helps with testing
    gradient(repmat(logical(obstacle), [1 1 2])) = 0; 
end
