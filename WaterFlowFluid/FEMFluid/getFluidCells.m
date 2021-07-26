function [ fluidCells ] = getFluidCells( distanceField )
%GETFLUIDCELLS Convert a distance field to a boolean matrix. Fluid cells 
%       must contain at least one fluid node.
%   return fluidCells: boolean matrix(yCells, xNodes)
%   distanceField: matrix(yCells, xNodes)

    xCells = size(distanceField, 1) - 1;
    yCells = size(distanceField, 2) - 1;
    
    % surrounded by 4 nodes
    topLeftFluidNodes = distanceFieldToFluidNodes(...
        distanceField(1 : xCells, 1 : yCells));
    bottomLeftFluidNodes = distanceFieldToFluidNodes(...
        distanceField(1 : xCells, 2 : yCells+1));
    topRightFluidNodes = distanceFieldToFluidNodes(...
        distanceField(2 : xCells+1, 1 : yCells));
    bottomRightFluidNodes = distanceFieldToFluidNodes(...
        distanceField(2 : xCells+1, 2 : yCells+1));
    
    % mark water filled cells
    % treat mixed cells like completely water filled
    fluidCells = topLeftFluidNodes | bottomLeftFluidNodes | ...
                   topRightFluidNodes | bottomRightFluidNodes;
            
end

