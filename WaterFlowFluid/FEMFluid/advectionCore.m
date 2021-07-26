function [ advectedSource ] = advectionCore( ...
    targetCoords, source, interpolationAlgorithm )
%ADVECTIONCORE This function advects a source value field to given
%       target coordinates. This function takes care of the bounds by
%       keeping the target coordinates within the value field size.
%   return advectedSource: matrix(yNodes, xNodes, [dim]), 
%       advected state of the source matrix, velocity or distance field
%   targetCoords: matrix(yNodes, xNodes)
%   source: matrix(yNodes, xNodes, [dim]), velocity or distance field
%   interpolationAlgorithm: string, one of these 'spline', 'cubic',
%       'linear', 'nearest'. Linear causes a major loss of fluid.
%       Spline minimizes the loss.
    
    mSize = size(source, 1);
    nSize = size(source, 2);
    
    % add padding, -1 because we don't want fluid to come out of nowhere
    source = padarray(source, [1 1], -1);
    targetCoords = targetCoords + 1;
    mPaddedSize = size(source, 1);
    nPaddedSize = size(source, 2);
    targetCoords = padarray(targetCoords, [1 1], 'replicate');
    
    % keep coordinates in valid range
    targetCoords(targetCoords < 1) = 1;
    xCoords = targetCoords(:, :, 1);
    yCoords = targetCoords(:, :, 2);
    xCoords(xCoords > nPaddedSize) = nPaddedSize;
    yCoords(yCoords > mPaddedSize) = mPaddedSize;
    targetCoords(:, :, 1) = xCoords;
    targetCoords(:, :, 2) = yCoords;
    
    % interpolate only 1 dim for distance field, or 2 independent dims for
    % velocities
    dimsToInterpolate = size(source, 3);
    advectedSource = zeros(size(source));
    for i = 1 : dimsToInterpolate
        % bilinear interpolation around origin of velocity
        advectedSource(:,:,i) = interp2(...
            source(:,:,i), targetCoords(:,:,1), targetCoords(:,:,2), ...
            interpolationAlgorithm);
    end

    advectedSource = advectedSource(2 : mSize + 1, 2 : nSize + 1, :);
    
end

