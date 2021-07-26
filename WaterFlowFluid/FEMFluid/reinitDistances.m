function [ newDistanceField ] = reinitDistances( distanceField, obstacle )
%REINITDISTANCES Re-initialize the distance field. Nodes at the interface 
%       keep their distance. Without this function, the interpolation in 
%       the advection would cause the distances to blend.
%   return newDistanceField: matrix(yNodes, xNodes)
%   distanceField: matrix(yNodes, xNodes)
    
    % get fluid cells and interface nodes
    interfaceNodes = getInterfaceNodes(distanceField, obstacle, false);
    
    % get water filled areas
    fluidNodes = distanceFieldToFluidNodes(distanceField);
    
    % assemble new distance field
    newDistanceField = fluidNodesToDistanceField(fluidNodes);
    
    % remove infinite values (created by bwdist if no water exists)
    newDistanceField(isinf(newDistanceField)) = -1;
    
    % keep distance at the interface
    newDistanceField(interfaceNodes) = distanceField(interfaceNodes);

end

