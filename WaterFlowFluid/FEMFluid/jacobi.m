function [ newPressure ] = jacobi( pressure, divergence, obstacle )
%JACOBI Calculate the pressure
% based on: http://scrawkblog.com/2013/05/21/gpu-gems-to-unity-2d-fluid-simulation/
%   return newPressure: matrix(yNodes, xNodes)
%   pressure: matrix(yNodes, xNodes)
%   divergences: matrix(yNodes, xNodes)
%   obstacle: boolean matrix(yNodes, xNodes)

    alpha = -1;
    inverseBeta = 0.25;

    xSize = size(pressure, 1);
    ySize = size(pressure, 2);
    
    % make look-up at borders easier
    % same effect as clamped texture coordinates in opengl
    paddedPressure = padarray(pressure, [1 1], 'replicate');
    paddedObstacle = padarray(obstacle, [1 1], 'replicate');
    
    % Find neighboring pressure
    pressureOneUp = paddedPressure(2 : xSize + 1, 3 : ySize + 2);
    pressureOneDown = paddedPressure(2 : xSize + 1, 1 : ySize);
    pressureOneRight = paddedPressure(3 : xSize + 2, 2 : ySize + 1);
    pressureOneLeft = paddedPressure(1 : xSize, 2 : ySize + 1);
    
    % Find neighboring obstacles
    obstacleOneUp = paddedObstacle(2 : xSize + 1, 3 : ySize + 2);
    obstacleOneDown = paddedObstacle(2 : xSize + 1, 1 : ySize);
    obstacleOneRight = paddedObstacle(3 : xSize + 2, 2 : ySize + 1);
    obstacleOneLeft = paddedObstacle(1 : xSize, 2 : ySize + 1);
    
    % Use center pressure for solid cells
    pressureOneUp(obstacleOneUp == true) = pressure(obstacleOneUp == true);
    pressureOneDown(obstacleOneDown == true) = pressure(obstacleOneDown == true);
    pressureOneRight(obstacleOneRight == true) = pressure(obstacleOneRight == true);
    pressureOneLeft(obstacleOneLeft == true) = pressure(obstacleOneLeft == true);
    
    newPressure = (pressureOneLeft + pressureOneRight + ...
        pressureOneDown + pressureOneUp + alpha .* divergence) .* inverseBeta;
end

