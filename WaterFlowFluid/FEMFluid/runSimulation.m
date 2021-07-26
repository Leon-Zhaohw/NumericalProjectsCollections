function [divergence, fluidPixels] = runSimulation( ...
    simulationSettings, visualizationSettings )
%RUNSIMULATION Run the simulation and save frames as images
%   return divergence: vector(numFrames - 1) sum of divergence
%   return fluidPixels: vector(numFrames - 1) num fluid pixels * 1000
%   simulationSettings: struct with these fields:
%       xCells: grid cells in x direction
%       yCells: grid cells in y direction
%       cellSize: grid cell size in m
%       noSteps: number of simulation steps as integer
%       deltaTime: passed time per step in seconds
%       advectionAlgorithm: string, one of these 'semilagrange', 
%           'maccormack'
%       advectionInterpolation: string, one of these 'spline', 'cubic', 
%           'linear', 'nearest'
%       dissipation: friction within fluid
%       scenario: string as defined in loadScenario
%   visualizationSettings: struct with these fields:
%       outputFolder: folder where the resulting video should be saved
%       gridCellSize: inflate cells to pixels
%       isobarWidth: width of the isobars (areas with similar distance)
%       drawVelocity: bool, display velocity field
%       distanceRefreshInterval: integer, distance field is reinitialized 
%           after n frames

    display('Loading scenario...');
    
    [ distanceField, velocityField, pressureField, ...
        obstacle, forceField] = loadScenario(simulationSettings);
    
    display('Simulating...');
    
    % select advection algorithm
    if strcmp(simulationSettings.advectionAlgorithm, 'maccormack')
        % add only constant parameters here as they won't change
        advect = @(velField, toBeAdvected) advectMacCormack( ...
            simulationSettings, velField, toBeAdvected);
    elseif strcmp(simulationSettings.advectionAlgorithm, 'semilagrange')
        advect = @(velField, toBeAdvected) advectSemiLagrange( ...
            simulationSettings, velField, toBeAdvected);
    else
        display('Wrong advection algorithm!');
        return;
    end
    
    % save start condition
    frame = calcFrame(visualizationSettings, distanceField, ...
        obstacle, velocityField, pressureField);
    saveImage(frame, 0, visualizationSettings.outputFolder);
        
% for debugging:
% quiver(flip(velocityField(:,:,1),1),-flip(velocityField(:,:,2),1))
% imshow(distanceField>=0)
% HeatMap(pressureField)
    
    % quality measures
    divergence = zeros(1, simulationSettings.numFrames - 1);
    fluidPixels = zeros(1, simulationSettings.numFrames - 1);
        
    % delete velocities in obstacles (set to obstacle velocity)
    velocityField(repmat(obstacle, [1 1 2])) = 0;
        
    for i = 1 : simulationSettings.numFrames - 1
        %% simulation
        
        % apply forces
        velocityField = applyForces(velocityField, forceField, simulationSettings.deltaTime);
        
        % delete velocities in obstacles (set to obstacle velocity)
        velocityField(repmat(obstacle, [1 1 2])) = 0;
        
        % apply pressure correction
        pressureField = solvePressure(simulationSettings, velocityField, distanceField, obstacle);
%         pressureField = calcDivergence(simulationSettings, velocityField, distanceField, obstacle);
        
        % pressure projection to keep the water volume (as much as possible)
        velocityField = subtractPressureGradient(velocityField, pressureField, obstacle);
        
        % correct velocities in air, necessary for advection
        velocityField = setAirVelocity(velocityField, distanceField, obstacle);
        
        % Advect density against velocity
        distanceField = advect(velocityField, distanceField);
        % Advect velocity against itself
        velocityField = advect(velocityField, velocityField);
        
        %% visualization
        % re-init distance field after n frames
        if visualizationSettings.distanceRefreshInterval > 0 && ...
                mod(i, visualizationSettings.distanceRefreshInterval) == 0
            distanceField = reinitDistances(distanceField, obstacle);
        end

        [frame, fluidPixels(i)] = calcFrame(visualizationSettings, distanceField, ...
            obstacle, velocityField, pressureField);
        saveImage(frame, i, visualizationSettings.outputFolder);
        
        % store overall divergence (should be 0 after perfect pressure solving)
        waterBody = distanceFieldToFluidNodes(distanceField);
        frameDivergence=calcDivergence( velocityField, obstacle );
        frameDivergence(obstacle)=0; % care only for fluid
        frameDivergence(~waterBody)=0;
        divergence(i) = sum(sum(frameDivergence));
    
        disp(strcat('Frame:', num2str(i)));
    end
    
    display('Simulation finished!');
end

