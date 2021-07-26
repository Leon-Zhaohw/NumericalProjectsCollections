function [ advectedSource ] = advectSemiLagrange( ...
    simulationSettings, velocity, source )
%ADVECTSEMILAGRANGE Apply Semi-Lagrangian advection on source data with 
%       the velocity field
%   return advectedSource: matrix(yNodes, xNodes, [dim]), 
%       advected state of the source matrix, velocity or distance field
%   simulationSettings: struct with these fields:
%       cellSize: grid cell size in m
%       deltaTime: passed time per step in seconds
%       dissipation: friction within fluid
%       advectionInterpolation: string, one of these 'spline', 'cubic', 
%           'linear', 'nearest'. Linear causes a major loss of fluid.
%           Spline minimizes the loss.
%   velocity: matrix(yNodes, xNodes, dim) in m/s
%   source: matrix(yNodes, xNodes, [dim]), velocity or distance field
    
    mSize = size(velocity, 1);
    nSize = size(velocity, 2);

    % create matrix, which contains coordinates from 1 to width / height
    coords = zeros(mSize, nSize, 2);
    [coords(:, :, 1), coords(:, :, 2)] = meshgrid(1 : nSize, 1 : mSize);
    
    % track coordinates back to their origin in the last step
    velocityInCellsPerSec = velocity / simulationSettings.cellSize;
    lookUpcoords = coords - velocityInCellsPerSec * simulationSettings.deltaTime;
    
    interpolatedSource = advectionCore(...
        lookUpcoords, source, simulationSettings.advectionInterpolation );
    
    % add dissipation (friction)
    advectedSource = simulationSettings.dissipation .* interpolatedSource;

end

