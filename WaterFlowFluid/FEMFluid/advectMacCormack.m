function [ advectedSource ] = advectMacCormack( ...
    simulationSettings, velocity, source )
%ADVECTMacCormack Apply advection based on MacCormack's Scheme on source 
%       data with the velocity field. More in Selle, Andrew, et al. 
%       "An unconditionally stable MacCormack method." Journal of 
%       Scientific Computing 35.2-3 (2008): 350-371,
%       http://physbam.stanford.edu/~fedkiw/papers/stanford2006-09.pdf
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
    
    % Calculation steps:
    % (1) normal semi-lagrangian advection from current state
    % (2) backward semi-lagrangian advection from (1)
    % (3) res = (1) + (currentState - (2)) / 2
    % (4) clamp with values around (1)

    % create matrix, which contains coordinates from 1 to width / height
    coords = zeros(mSize, nSize, 2);
    [coords(:, :, 1), coords(:, :, 2)] = meshgrid(1 : nSize, 1 : mSize);
    
    velocityInCellsPerSec = velocity / simulationSettings.cellSize;
    
    % (1) normal semi-lagrangian advection from current state
    forwardStepCoords = coords - (velocityInCellsPerSec * simulationSettings.deltaTime);
    
    % (2) normal semi-lagrangian advection from current state
    backwardStepCoords = forwardStepCoords + (velocityInCellsPerSec * simulationSettings.deltaTime);
    
    % (3) res = (1) + (currentState - (2)) / 2
    lookUpcoords = forwardStepCoords + ...
        (coords - backwardStepCoords) ./ 2;
    
    % interpolate only 1 dim for distance field, or 2 independent dims
    % for velocities
    finalSource = advectionCore(lookUpcoords, source, ...
        simulationSettings.advectionInterpolation );
    
    % interpolate for clamping
    forwardSource = advectionCore(forwardStepCoords, source, ...
        simulationSettings.advectionInterpolation );
    
    % (4) clamp with values around (1)
    % interpolate only 1 dim for distance field, or 2 independent dims for
    % velocities
    numSourceDims = size(source, 3);
    for i = 1 : numSourceDims
        % get maximum values around each cell, excluding center
        neighborhood = [1 1 1; 1 0 1; 1 1 1];
        maxValues = imdilate(forwardSource(:,:,i), neighborhood);
        minValues = imerode(forwardSource(:,:,i), neighborhood);
        
        % apply clamp
        finalSourceDim = finalSource(:,:,i);
        finalSourceDim(finalSourceDim > maxValues) = maxValues(finalSourceDim > maxValues);
        finalSourceDim(finalSourceDim < minValues) = minValues(finalSourceDim < minValues);
        finalSource(:,:,i) = finalSourceDim;
    end
    
    % add dissipation (friction)
    advectedSource = simulationSettings.dissipation .* finalSource;

end

