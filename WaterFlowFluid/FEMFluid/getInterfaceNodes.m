function [ interfaceNodes ] = getInterfaceNodes(distanceField, obstacle, onlyAir)
%GETINTERFACENODES Get nodes of cells in which the interface is. Interfaces
%       at obstacles don't count as interfaces.
%   return interfaceNodes: boolean matrix(yNodes, xNodes)
%   distanceField: matrix(yNodes, xNodes)
%   obstacle: boolean matrix(yNodes, xNodes)
%   onlyAir: boolean scalar

ySize = size(distanceField, 1);
xSize = size(distanceField, 2);

% in case someone forgot to make it logical
obstacle = logical(obstacle);

fluid = distanceFieldToFluidNodes(distanceField);
paddedFluid = padarray(fluid, [1 1]);
paddedObstacle = padarray(obstacle, [1 1]);

%% Find neighboring fluid nodes
% neighborhood of 4
fluidLeft = logical(paddedFluid(2 : ySize + 1, 3 : xSize + 2));
fluidRight = logical(paddedFluid(2 : ySize + 1, 1 : xSize));
fluidUp = logical(paddedFluid(3 : ySize + 2, 2 : xSize + 1));
fluidDown = logical(paddedFluid(1 : ySize, 2 : xSize + 1));
% neighborhood of 8
fluidDownLeft = logical(paddedFluid(1 : ySize, 3 : xSize + 2));
fluidUpLeft = logical(paddedFluid(3 : ySize + 2, 3 : xSize + 2));
fluidDownRight = logical(paddedFluid(1 : ySize, 1 : xSize));
fluidUpRight = logical(paddedFluid(3 : ySize + 2, 1 : xSize));

%% Find neighboring obstacles
% neighborhood of 4
obstacleLeft = logical(paddedObstacle(2 : ySize + 1, 3 : xSize + 2));
obstacleRight = logical(paddedObstacle(2 : ySize + 1, 1 : xSize));
obstacleUp = logical(paddedObstacle(3 : ySize + 2, 2 : xSize + 1));
obstacleDown = logical(paddedObstacle(1 : ySize, 2 : xSize + 1));
% neighborhood of 8
obstacleDownLeft = logical(paddedObstacle(1 : ySize, 3 : xSize + 2));
obstacleUpLeft = logical(paddedObstacle(3 : ySize + 2, 3 : xSize + 2));
obstacleDownRight = logical(paddedObstacle(1 : ySize, 1 : xSize));
obstacleUpRight = logical(paddedObstacle(3 : ySize + 2, 1 : xSize));


%% a node is interface node if:
% - it is not in an obstacle
% - at least one neighbor's fluid state is different than the middle
% - that neighbor node is not in an obstacle
interfaceNodes = false(size(obstacle));

% neighborhood of 4
%     interfaceNodes = interfaceNodes | ...
%         (fluid ~= fluidUp & ~obstacleUp) |  ...
%         (fluid ~= fluidDown & ~obstacleDown) | ...
%         (fluid ~= fluidRight & ~obstacleRight) | ...
%         (fluid ~= fluidLeft & ~obstacleLeft);
interfaceNodes = interfaceNodes | (fluid ~= fluidUp & ~obstacleUp);
interfaceNodes = interfaceNodes |(fluid ~= fluidDown & ~obstacleDown);
interfaceNodes = interfaceNodes | (fluid ~= fluidRight & ~obstacleRight);
interfaceNodes = interfaceNodes | (fluid ~= fluidLeft & ~obstacleLeft);
% neighborhood of 8
interfaceNodes = interfaceNodes | ...
    (fluid ~= fluidUpLeft & ~obstacleUpLeft) | ...
    (fluid ~= fluidUpRight & ~obstacleUpRight) | ...
    (fluid ~= fluidDownLeft & ~obstacleDownLeft) | ...
    (fluid ~= fluidDownRight & ~obstacleDownRight);

% center must not be in an obstacle
interfaceNodes(obstacle) = false;

% experiment: take only air nodes for interface
if onlyAir
    interfaceNodes(fluid) = false;
end

end