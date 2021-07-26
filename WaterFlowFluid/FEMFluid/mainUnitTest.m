function mainUnitTest()

clear global; % no clear all or breakpoints are gone
close all;
clc;

% break on error
dbstop if error

numErrors = 0;
disp(strcat('Start unit test at: ', datestr(now)));

error = testSubtractPressureGradient();
if error ~= 0
    numErrors = numErrors + 1;
end

error = testExpansion();
if error ~= 0
    numErrors = numErrors + 1;
end

error = testBoundaries();
if error ~= 0
    numErrors = numErrors + 1;
end

disp(strcat('Finish unit test at: ', datestr(now)));
disp(strcat('Number of failed tests: ', num2str(numErrors)));
end

function [ error ] = testSubtractPressureGradient()

    obstacle = ...
    [
        1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
        1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1;
        1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
        1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
    ];

    pressureField = -1 .* ...
    [
        0    0   0   0   0   0   0   0   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0    0   0   0   0   0   0   0   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  1.4   0   0   0   0   0   0   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  1.9   0   0   0   0   0   0   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.1 1.5 0.7   0   0   0   0   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.2 1.7 1.1 0.8 0.5   0   0   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.2 1.8 1.4 1.1 0.8 0.5   0   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.2 1.9 1.6 1.3 1.0 0.7 0.4   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.2 2.0 1.7 1.5 1.2 0.9 0.5 0.3   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.3 2.0 1.8 1.6 1.3 1.0 0.7 0.4 0.2   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.3 2.1 1.8 1.7 1.4 1.1 0.8 0.5 0.3   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.3 2.1 1.8 1.7 1.4 1.1 0.9 0.6 0.4 0.2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.3 2.1 1.9 1.7 1.4 1.2 0.9 0.6 0.4 0.2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.3 2.1 1.8 1.7 1.4 1.1 0.9 0.6 0.4 0.2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.3 2.1 1.8 1.7 1.4 1.1 0.8 0.5 0.3   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.3 2.0 1.8 1.6 1.3 1.0 0.7 0.4 0.2   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.2 2.0 1.7 1.5 1.2 0.9 0.5 0.3   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.2 1.9 1.6 1.3 1.0 0.7 0.4   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.2 1.8 1.4 1.1 0.8 0.5   0   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.2 1.7 1.1 0.8 0.5   0   0   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  2.1 1.5 0.7   0   0   0   0   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  1.9   0   0   0   0   0   0   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  1.4   0   0   0   0   0   0   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0    0   0   0   0   0   0   0   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0    0   0   0   0   0   0   0   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0    0   0   0   0   0   0   0   0   0   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    ]';

    expectedGradient = cat(3, ...
    [
        0  0     0     0     0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0     0     0     0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0     0     0     0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0     0     0     0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.80 -0.75 -0.35     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.60 -0.45 -0.30 -0.40 -0.25     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.40 -0.35 -0.30 -0.30 -0.40 -0.25     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.30 -0.30 -0.30 -0.30 -0.30 -0.35 -0.20     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.30 -0.25 -0.25 -0.30 -0.35 -0.30 -0.25 -0.15     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.20 -0.20 -0.25 -0.30 -0.30 -0.30 -0.25 -0.20 -0.10    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.30 -0.20 -0.20 -0.30 -0.30 -0.30 -0.25 -0.25 -0.15    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.30 -0.20 -0.20 -0.30 -0.25 -0.25 -0.25 -0.20 -0.20 -0.1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.20 -0.20 -0.25 -0.25 -0.25 -0.30 -0.25 -0.20 -0.20 -0.1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.30 -0.20 -0.20 -0.30 -0.25 -0.25 -0.25 -0.20 -0.20 -0.1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.30 -0.20 -0.20 -0.30 -0.30 -0.30 -0.25 -0.25 -0.15    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.20 -0.20 -0.25 -0.30 -0.30 -0.30 -0.25 -0.20 -0.10    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.30 -0.25 -0.25 -0.30 -0.35 -0.30 -0.25 -0.15     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.30 -0.30 -0.30 -0.30 -0.30 -0.35 -0.20     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.40 -0.35 -0.30 -0.30 -0.40 -0.25     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.60 -0.45 -0.30 -0.40 -0.25     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0 -0.80 -0.75 -0.35     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0     0    0      0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0     0    0      0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0     0    0      0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0     0    0      0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0  0     0    0      0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    ]', ...
    [
        0 0     0     0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0     0     0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0     0     0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0  0.75  0.35     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0  0.85  0.55  0.40  0.25     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0  0.15  0.35  0.55  0.40  0.25     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0  0.10  0.25  0.25  0.25  0.35  0.20     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0  0.10  0.15  0.20  0.20  0.20  0.25  0.15     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0  0.05  0.10  0.15  0.15  0.15  0.15  0.20  0.10    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0  0.05  0.05  0.10  0.10  0.10  0.15  0.10  0.15    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0  0.05     0  0.05  0.05  0.05  0.10  0.10  0.10  0.1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0     0  0.05     0     0  0.05  0.05  0.05  0.05  0.1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0     0     0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0     0 -0.05     0     0 -0.05 -0.05 -0.05 -0.05 -0.1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -0.05     0 -0.05 -0.05 -0.05 -0.10 -0.10 -0.10 -0.1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -0.05 -0.05 -0.10 -0.10 -0.10 -0.15 -0.10 -0.15    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -0.05 -0.10 -0.15 -0.15 -0.15 -0.15 -0.20 -0.10    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -0.10 -0.15 -0.20 -0.20 -0.20 -0.25 -0.15     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -0.10 -0.25 -0.25 -0.25 -0.35 -0.20     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -0.15 -0.35 -0.55 -0.40 -0.25     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -0.85 -0.55 -0.40 -0.25     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -0.75 -0.35     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0     0     0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0     0     0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0     0     0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0     0     0     0     0     0     0     0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    ]');

    myEps = 0.001;

    gradient = calcGradient(pressureField, obstacle);
    
    errors = abs(expectedGradient - gradient) > myEps;
    if sum(sum(sum(errors))) == 0 
        disp('Successful: testSubtractPressureGradient');
        error = 0;
    else
        disp('Failed: testSubtractPressureGradient');
        disp('Wrong elements:');
        disp(errors);
        error = 1;
    end
    
end



function [ error ] = testExpansion()

    %% initialize
    xCells = 18;
    yCells = 9;

    domainSize = 1; %m
    simulationTime = 3; %s
    deltaTime = 1 / 60;

    simulationSettings = struct(...
        'xCells', xCells, ...
        'yCells', yCells, ...
        'cellSize', domainSize / xCells, ...
        'numFrames', simulationTime / deltaTime, ...
        'deltaTime', deltaTime, ...
        'advectionAlgorithm', 'maccormack', ...
        'dissipation', 1, ...
        'scenario', 'expandingBlob' ...
    );
    
    [ distanceField, velocityField, ~, ...
        obstacle, ~] = loadScenario(simulationSettings);
    
    expectedVelocity = zeros(yCells + 1, xCells + 1, 2);

    myEps = 0.001;
    
    %% calc result
    pressureField = solvePressure(simulationSettings, velocityField, distanceField, obstacle);
    velocityField = subtractPressureGradient(velocityField, pressureField, obstacle);
    velocityField = enforceBoundaryCondition(velocityField, obstacle);
    
    %% compare result
    errors = abs(expectedVelocity - velocityField) > myEps;
    if sum(sum(sum(errors))) == 0 
        disp('Successful: testExpansion');
        error = 0;
    else
        disp('Failed: testExpansion');
        disp('Wrong elements:');
        disp(errors);
        error = 1;
    end
    
end

function [ error ] = testBoundaries()

    %% initialize
    xCells = 9;
    yCells = 9;

    domainSize = 1; %m
    simulationTime = 3; %s
    deltaTime = 1 / 60;

    simulationSettings = struct(...
        'xCells', xCells, ...
        'yCells', yCells, ...
        'cellSize', domainSize / xCells, ...
        'numFrames', simulationTime / deltaTime, ...
        'deltaTime', deltaTime, ...
        'advectionAlgorithm', 'maccormack', ...
        'dissipation', 1, ...
        'scenario', 'restingBlob' ...
    );
    
    [ distanceField, velocityField, ~, ...
        obstacle, ~] = loadScenario(simulationSettings);
    
    expectedVelocity = zeros(yCells + 1, xCells + 1, 2);

    myEps = 0.001;
    
    %% calc result
%     velocityField = applyForces(velocityField, forceField, simulationSettings.deltaTime);
    pressureField = solvePressure(simulationSettings, velocityField, distanceField, obstacle);
    velocityField = subtractPressureGradient(velocityField, pressureField, obstacle);
    velocityField = enforceBoundaryCondition(velocityField, obstacle);
    
    %% compare result
    errors = abs(expectedVelocity - velocityField) > myEps;
    if sum(sum(sum(errors))) == 0 
        disp('Successful: testBoundaries');
        error = 0;
    else
        disp('Failed: testBoundaries');
        disp('Wrong elements:');
        disp(errors);
        error = 1;
    end
    
end