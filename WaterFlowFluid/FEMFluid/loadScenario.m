function [ distanceField, velocityField, pressureField, ...
    obstacle, forceField ] = ...
    loadScenario( simulationSettings )
%LOADSCENARIO Create blobs and obstacles as defined in a scenario
%   return distanceField: matrix(yNodes, xNodes), distances to interface
%       velocityField: matrix(yNodes, xNodes, dim), velocity of
%           fluid in m/s
%       pressureField: matrix(yNodes, xNodes), pressure at each node
%       obstacle: bool matrix(yNodes, xNodes), node is in obstacle 
%           (no fluid or air)
%       forceField: matrix(yNodes, xNodes, dimension), force applied to
%           fluid in m/s^2
%   simulationSettings: struct with these fields:
%       xCells: grid cells in x direction
%       yCells: grid cells in y direction
%       cellSize: grid cell size in m
%       scenario: string, one of the construction functions below
%       deltaTime: passed time per step in seconds
    
    % never initial pressure
    xNodes = simulationSettings.xCells + 1;
    yNodes = simulationSettings.yCells + 1;
    pressureField = zeros(yNodes, xNodes);
    
    try
        functionHandle = str2func(simulationSettings.scenario);
        [ fluidNodes, velocityField, obstacle, forceField] = ...
            functionHandle(simulationSettings);
    catch ME
        switch ME.identifier
            case 'MATLAB:UndefinedFunction'
                warning(strcat('Invalid scnenario "', simulationSettings.scenario, ...
                    '" choosen. Use the name of a function in loadScenario.m.', ...
                    'Loading scenario "twoBlobs" instead'));
            otherwise
                warning(strcat('Invalid scnenario "', simulationSettings.scenario, ...
                    '" choosen. Use the name of a function in loadScenario.m.', ...
                    'Loading scenario "twoBlobs" instead. ', ...
                    'This message can mean that the chosen scenario has an error. ', ...
                    'Error: ', ME.identifier, ' Message: ', ME.message));
        end
        
        [ fluidNodes, velocityField, obstacle, forceField] = ...
            twoCollidingBlobs(simulationSettings);
    end
    
    % cells to nodes
    obstacle = logical(cellsToNodes(obstacle));
    fluidNodes = cellsToNodes(fluidNodes);
    
    % calculate distanceField from waterbody
    distanceField = fluidNodesToDistanceField(fluidNodes);

end


function waterBody = makeRoundBlob( xCells, yCells, ...
    blobCenterXPercentage, blobCenterYPercentage, radiusPercentage )
    
    % round blob in middle
    [x, y] = meshgrid(1 : xCells, 1 : yCells);
    x = double(x);
    y = double(y);
    x = x - blobCenterXPercentage * xCells;
    y = y - blobCenterYPercentage * yCells;
    distToCenter = sqrt(x .* x + y .* y);
    waterBody = false(yCells, xCells);
    waterBody(distToCenter < radiusPercentage * min(xCells, yCells)) = true;
    
end

function [ waterBody, velocityField, obstacle, forceField] = ...
    advectionTest( simulationSettings )

    xCells = simulationSettings.xCells;
    yCells = simulationSettings.yCells;
    xNodes = simulationSettings.xCells + 1;
    yNodes = simulationSettings.yCells + 1;
    
    % move everything down
    velocityField = zeros(yNodes, xNodes, 2);
    velocityField(:, :, 2) = 0.4;
    
    % no forces
    forceField = zeros(yNodes, xNodes, 2);
    
    % one blob at top
    leftPercentage = 0.4;
    rightPercentage = 0.5;
    topPercentage = 0.2;
    bottomPercentage = 0.1;
    waterBody = false(yCells, xCells);
    waterBody(floor(yCells * bottomPercentage) : ...
        floor(yCells * topPercentage), ...
        floor(xCells * leftPercentage) : ...
        floor(xCells * rightPercentage)) = true;
    
    % scale as obstacles, every 10th cell
    obstacle = false(yCells, xCells);
    obstacle(1 : 10 : yCells, xCells) = true;
    % boundary as pipe with bottom
    obstacle(:, floor(xCells * leftPercentage) - 2) = true;
    obstacle(:, floor(xCells * rightPercentage) + 2) = true;
    obstacle(yCells, floor(xCells * leftPercentage) - 2 : ...
        floor(xCells * rightPercentage) + 2) = true;
end

function [ waterBody, velocityField, obstacle, forceField] = ...
    twoCollidingBlobs( simulationSettings )

    xCells = simulationSettings.xCells;
    yCells = simulationSettings.yCells;
    xNodes = simulationSettings.xCells + 1;
    yNodes = simulationSettings.yCells + 1;
    
    % move top right against left bottom
    velocityField = zeros(yNodes, xNodes, 2);
    
    velocityField(uint32(1 : yNodes * 0.5), uint32(1 : xNodes * 0.5), 1) = 0.2;
    velocityField(uint32(yNodes * 0.5 : yNodes), 1 : xNodes, 1) = -0.2;
    velocityField(uint32(1 : yNodes * 0.5), uint32(1 : xNodes * 0.5), 2) = 0.2;
    velocityField(uint32(yNodes * 0.5 : yNodes), 1 : xNodes, 2) = -0.2;
    
    % no forces
    forceField = zeros(yNodes, xNodes, 2);
    
    % separated blobs top right and bottom left
    waterBody = false(yCells, xCells);
    waterBody(uint32(yCells * 0.1 : yCells * 0.3), ...
        uint32(xNodes * 0.1 : xNodes * 0.3)) = true;
    waterBody(uint32(yCells * 0.55 : yCells * 0.75), ...
        uint32(xNodes * 0.55 : xNodes * 0.75)) = true;
    
    % obstacles: matrix(step, xPos, yPos)
    obstacleBorderWidth = 1;
    obstacle = false(yCells, xCells);
    
    % border around the fluid
    obstacle(1 : obstacleBorderWidth, :) = true;
    obstacle(yCells - obstacleBorderWidth + 1 : yCells, :) = true;
    obstacle(:, 1 : obstacleBorderWidth) = true;
    obstacle(:, xCells - obstacleBorderWidth + 1 : xCells) = true;
end

function [ waterBody, velocityField, obstacle, forceField] = ...
    stretchingBlob( simulationSettings )

    xCells = simulationSettings.xCells;
    yCells = simulationSettings.yCells;
    xNodes = simulationSettings.xCells + 1;
    yNodes = simulationSettings.yCells + 1;
    
    % linear gradient
    slope = 1 / double(yNodes);
    gradientLine = (-0.5 + 0.5 * slope : slope : 0.5 - 0.5 * slope)';
    velocityField = zeros(yNodes, xNodes, 2);
    velocityField(1 : yNodes, 1 : xNodes, 2) = repmat(gradientLine, [1 xNodes]);
    
    % no forces
    forceField = zeros(yNodes, xNodes, 2);
    
    % round blob in middle
    waterBody = makeRoundBlob( xCells, yCells, 0.5, 0.5, 0.2 );
    
    % obstacles: matrix(step, xPos, yPos)
    obstacleBorderWidth = 1;
    obstacle = false(yCells, xCells);
    
    % border around the fluid
    obstacle(1 : obstacleBorderWidth, :) = true;
    obstacle(yCells - obstacleBorderWidth + 1 : yCells, :) = true;
    obstacle(:, 1 : obstacleBorderWidth) = true;
    obstacle(:, xCells - obstacleBorderWidth + 1 : xCells) = true;
    
end

function [ waterBody, velocityField, obstacle, forceField] = ...
    expandingBlob( simulationSettings )

    xCells = simulationSettings.xCells;
    yCells = simulationSettings.yCells;
    xNodes = simulationSettings.xCells + 1;
    yNodes = simulationSettings.yCells + 1;
    
    % velocity radial from center
    [x, y] = meshgrid(1 : xNodes, 1 : yNodes);
    x = double(x);
    y = double(y);
    x = x - 0.5 * xNodes;
    y = y - 0.5 * yNodes;
    velocityField = cat(3, x, y);
    velocityField = velocityField ./ simulationSettings.yCells;
    
    % no forces
    forceField = zeros(yNodes, xNodes, 2);
    
    % round blob in middle
    waterBody = makeRoundBlob( xCells, yCells, 0.5, 0.5, 0.3 );
    
    % no obstacle
    obstacle = zeros(yCells, xCells);
end

function [ waterBody, velocityField, obstacle, forceField] = ...
    singleCollidingBlob( simulationSettings )

    xCells = simulationSettings.xCells;
    yCells = simulationSettings.yCells;
    xNodes = simulationSettings.xCells + 1;
    yNodes = simulationSettings.yCells + 1;
    
    % move top against bottom
    velocityField = zeros(yNodes, xNodes, 2);
    velocityField(1 : yNodes, uint32(xNodes * 0.5 : xNodes), 1) = -0.2;
    velocityField(1 : yNodes, uint32(1 : xNodes * 0.5), 1) = 0.1;
    
    % no forces
    forceField = zeros(yNodes, xNodes, 2);
    
    % big blob
    waterBody = false(yCells, xCells);
    waterBody(uint32(yCells * 0.25 : yCells * 0.75), ...
        uint32(xCells * 0.25 : xCells * 0.75), :) = true;
    
    % obstacles: matrix(step, xPos, yPos)
    obstacleBorderWidth = 1;
    obstacle = false(yCells, xCells);
    
    % border around the fluid
    obstacle(1 : obstacleBorderWidth, :) = true;
    obstacle(yCells - obstacleBorderWidth + 1 : yCells, :) = true;
    obstacle(:, 1 : obstacleBorderWidth) = true;
    obstacle(:, xCells - obstacleBorderWidth + 1 : xCells) = true;
    
end

function [ waterBody, velocityField, obstacle, forceField] = ...
    singleBlocker( simulationSettings )

    xCells = simulationSettings.xCells;
    yCells = simulationSettings.yCells;
    xNodes = simulationSettings.xCells + 1;
    yNodes = simulationSettings.yCells + 1;
    
    % move everything down
%     velocityField = zeros(yNodes, xNodes, 2);
    velocityField(1 : yNodes, 1 : xNodes, 2) = 0.5;
    
    % no forces
    forceField = zeros(yNodes, xNodes, 2);
%     % apply gravity
%     forceField = zeros(yNodes, xNodes, 2);
    forceField(:, :, 2) = 9.81;
    
    % blob above box
    waterBody = makeRoundBlob( xCells, yCells, 0.5, 0.3, 0.25 );
    
    % obstacles: matrix(step, xPos, yPos)
    obstacleBorderWidth = 1;
    obstacle = false(yCells, xCells);
    % border around the fluid
%     obstacle(1 : obstacleBorderWidth, :) = true;
    obstacle(yCells - obstacleBorderWidth + 1 : yCells, :) = true;
    obstacle(:, 1 : obstacleBorderWidth) = true;
    obstacle(:, xCells - obstacleBorderWidth + 1 : xCells) = true;
    % box in middle
    obstacle(uint32(yCells * 0.65 : yCells * 0.85), ...
        uint32(xCells * 0.4 : xCells * 0.6)) = true;
    
end

function [ waterBody, velocityField, obstacle, forceField] = ...
    singleBlobInTrap( simulationSettings )

    xCells = simulationSettings.xCells;
    yCells = simulationSettings.yCells;
    xNodes = simulationSettings.xCells + 1;
    yNodes = simulationSettings.yCells + 1;
    
    % move everything down and left
    velocityField = zeros(yNodes, xNodes, 2);
    velocityField(1 : yNodes, 1 : xNodes, 1) = 0.5;
    velocityField(1 : yNodes, 1 : xNodes, 2) = 0.5;
    
    % no forces
    forceField = zeros(yNodes, xNodes, 2);
    
    % blob above trap
    waterBody = false(yCells, xCells);
    waterBody(uint32(yCells * 0.15 : yCells * 0.45), ...
        uint32(xCells * 0.25 : xCells * 0.55), :) = true;
    
    % obstacles: matrix(step, xPos, yPos)
    obstacleBorderWidth = 1;
    obstacle = false(yCells, xCells);
    
    % border around the fluid
    obstacle(1 : obstacleBorderWidth, :) = true;
    obstacle(yCells - obstacleBorderWidth + 1 : yCells, :) = true;
    obstacle(:, 1 : obstacleBorderWidth) = true;
    obstacle(:, xCells - obstacleBorderWidth + 1 : xCells) = true;
    
end

function [ waterBody, velocityField, obstacle, forceField] = ...
    singleBlobInTrapWithBlocker( simulationSettings )

    xCells = simulationSettings.xCells;
    yCells = simulationSettings.yCells;
    xNodes = simulationSettings.xCells + 1;
    yNodes = simulationSettings.yCells + 1;
    
    % move everything down and left
    velocityField = zeros(yNodes, xNodes, 2);
    velocityField(1 : yNodes, 1 : xNodes, 1) = -0.75;
%     velocityField(1 : yNodes, 1 : xNodes, 2) = 0.75;
    
%     % no forces
%     forceField = zeros(yNodes, xNodes, 2);
    % apply gravity
    forceField = zeros(yNodes, xNodes, 2);
    forceField(:, :, 2) = 9.81;
    
    % blob above trap
    waterBody = makeRoundBlob( xCells, yCells, 0.7, 0.3, 0.2 );
    
    % obstacles: matrix(step, xPos, yPos)
    obstacleBorderWidth = 1;
    obstacle = false(yCells, xCells);
    
    % border around the fluid
    obstacle(1 : obstacleBorderWidth, :) = true;
    obstacle(yCells - obstacleBorderWidth + 1 : yCells, :) = true;
    obstacle(:, 1 : obstacleBorderWidth) = true;
    obstacle(:, xCells - obstacleBorderWidth + 1 : xCells) = true;
    
    % blocker
    obstacle(uint32(yCells * 0.55 : yCells * 0.65), ...
        uint32(xCells * 0.35 : xCells * 0.45)) = true;    % bottom
    
end

function [ waterBody, velocityField, obstacle, forceField] = ...
    restingBlob( simulationSettings )

    xCells = simulationSettings.xCells;
    yCells = simulationSettings.yCells;
    xNodes = simulationSettings.xCells + 1;
    yNodes = simulationSettings.yCells + 1;
    
    % move down
    velocityField = zeros(yNodes, xNodes, 2);
%     velocityField(1 : yNodes, 1 : xNodes, 1) = -0.5;

%     % no force
%     forceField = zeros(yNodes, xNodes, 2);
    % apply gravity
    forceField = zeros(yNodes, xNodes, 2);
    forceField(:, :, 2) = 9.81;
    
    % obstacles: matrix(step, xPos, yPos)
    obstacleBorderWidth = 1;
    obstacle = false(yCells, xCells);
    
    % border around the fluid with padding
    padding = max(floor(yCells / 10), 1);
    obstacle(:, 1 + padding : obstacleBorderWidth + padding) = true;
    obstacle(:, xCells - obstacleBorderWidth + 1 - padding : xCells - padding) = true;
    obstacle(yCells - (1 + padding) : yCells - (obstacleBorderWidth + padding), :) = true;
    
    % blob in middle
    waterBody = false(yCells, xCells);
    waterBody(2 + obstacleBorderWidth + padding : yCells - padding - 3, ...
        2 + obstacleBorderWidth + padding : xCells - padding - 2) = true;
    
end

function [ waterBody, velocityField, obstacle, forceField] = ...
    fallingBlob( simulationSettings )

    xCells = simulationSettings.xCells;
    yCells = simulationSettings.yCells;
    xNodes = simulationSettings.xCells + 1;
    yNodes = simulationSettings.yCells + 1;
    
    % no initial movement
    velocityField = zeros(yNodes, xNodes, 2);
%     velocityField(1 : yNodes, 1 : xNodes, 2) = 1.0;
    % apply gravity
    forceField = zeros(yNodes, xNodes, 2);
    forceField(:, :, 2) = 9.81;
    
    % round blob in middle
    waterBody = makeRoundBlob( xCells, yCells, 0.5, 0.5, 0.25 );
    
    % obstacles: matrix(step, xPos, yPos)
    obstacleBorderWidth = 1;
    obstacle = false(yCells, xCells);
    
    % border around the fluid with padding
    padding = int32(yCells / 10);
    obstacle(:, 1 + padding : obstacleBorderWidth + padding) = true;
    obstacle(:, xCells - obstacleBorderWidth + 1 - padding : xCells - padding) = true;
    obstacle(yCells - (1 + padding) : yCells - (obstacleBorderWidth + padding), :) = true;
    
end

function [ waterBody, velocityField, obstacle, forceField] = ...
    manyBlockers ( simulationSettings )

    xCells = simulationSettings.xCells;
    yCells = simulationSettings.yCells;
    xNodes = simulationSettings.xCells + 1;
    yNodes = simulationSettings.yCells + 1;
    
    % no initial movement
    velocityField = zeros(yNodes, xNodes, 2);
    % apply gravity
    forceField = zeros(yNodes, xNodes, 2);
    forceField(:, :, 2) = 9.81;
    
    % round blob in middle
    waterAreaYSize = uint32(0.4 * yCells);
    waterBody = false(yCells, xCells);
    waterBody(2 : waterAreaYSize, 4 : xCells - 3) = true;
    
    % obstacles: matrix(xPos, yPos)
    obstacle = false(yCells, xCells);
    
    % many small interleaved obstacles
    smallPattern = [[1 1 1 0 1 1 1 0 1 1 1]; ...
                    [0 0 0 0 0 0 0 0 0 0 0]; ...
                    [1 0 1 1 1 0 1 1 1 0 1]; ...
                    [0 0 0 0 0 0 0 0 0 0 0]];
    smallPattern = double(repmat(smallPattern, [3 1]));
    [x, y] = meshgrid(1 : xCells, 1 : yCells);
    x = x / xCells * size(smallPattern, 2);
    y = y / yCells * size(smallPattern, 1);
    inflatedPattern = interp2(smallPattern, x, y, 'linear', 1);
    inflatedPattern(inflatedPattern < 0.65) = 0;
    inflatedPattern(inflatedPattern >= 0.65) = 1;
    obstacle(logical(inflatedPattern)) = true;
    
    % nothing in the fluid area
    obstacle(1 : waterAreaYSize + 2, :) = false;
    
    % border around the fluid with padding
    obstacle(:, 1) = true;
    obstacle(:, xCells) = true;
    obstacle(yCells, :) = true;
    
end

function [ waterBody, velocityField, obstacle, forceField] = ...
    funnel( simulationSettings )

    xCells = simulationSettings.xCells;
    yCells = simulationSettings.yCells;
    xNodes = simulationSettings.xCells + 1;
    yNodes = simulationSettings.yCells + 1;
    
    % no initial movement
    velocityField = zeros(yNodes, xNodes, 2);
    % apply gravity
    forceField = zeros(yNodes, xNodes, 2);
    forceField(:, :, 2) = 9.81;
    
    % round blob in middle
    waterBody = makeRoundBlob( xCells, yCells, 0.5, 0.3, 0.25 );
    
    % obstacles: matrix(xPos, yPos)
    obstacle = false(yCells, xCells);
    
    % funnel-like obstacle
    [x, y] = meshgrid(1 : xCells, 1 : yCells);
    x = double(x);
    y = double(y);
    x = x - 0.5 * xCells;
    y = y - 0.3 * yCells;
    distToCenter = abs(x) + y;
    obstacle((distToCenter > 0.7 * min(xCells, yCells)) & y >= 0) = true;
    
    % border around the fluid with padding
    padding = int32(yCells / 7.5);
    obstacle(:, 1 : padding) = true;
    obstacle(:, xCells - 1 - padding : xCells) = true;
    obstacle(yCells, :) = true;
    
end

function [ waterBody, velocityField, obstacle, forceField] = ...
    fallingCells( simulationSettings )

    xCells = simulationSettings.xCells;
    yCells = simulationSettings.yCells;
    xNodes = simulationSettings.xCells + 1;
    yNodes = simulationSettings.yCells + 1;
    
    % fluid size
    fluidSize = 0;
    
    % fluid position
    xFluid = 3;
    yFluid = yCells - 2;
    
    % move down 1 node per frame
    velocityField = zeros(yNodes, xNodes, 2);
    velocityField(1 : yNodes, 1 : xNodes, 2) = ...
        simulationSettings.cellSize / simulationSettings.deltaTime;

    % no force
    forceField = zeros(yNodes, xNodes, 2);
    
    % single cell in middle
    waterBody = false(yCells, xCells);
    waterBody(yFluid - fluidSize : yFluid, ...
        xFluid : xFluid + fluidSize) = true;
    
    % border as channel around the fluid
    obstacle = false(yCells, xCells);
    obstacle(:, xFluid - 2) = true;
    obstacle(:, xFluid + 2 + fluidSize) = true;
    obstacle(yCells, :) = true;
    
end
