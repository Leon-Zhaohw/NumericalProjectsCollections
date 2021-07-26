% Fluid Simulation with Finite Element Method
% As practical project for the visual computing master at the TU Wien
% 
% Written by Philipp Erler in 2016
% Supervised by Christian Hafner
% Feel free to ask, give feedback or just a note: ph.erler@gmx.net

% clear and close
clear global; % no clear all or breakpoints are gone
close all;
clc;

% break on error
dbstop if error

% profiling settings
enableProfiling = false;

% start new profiling session (clear old data)
if enableProfiling
    profile on
end

%'advectionAlgorithm': 'semilagrange', 'maccormack'
%'scenario': 'advectionTest', 'twoCollidingBlobs', 'stretchingBlob', ...
%   'expandingBlob', 'singleCollidingBlob', ...
%   'singleBlobInTrapWithBlocker', 'singleBlobInTrap', ...
%   'singleBlocker', 'restingBlob', 'fallingBlob', ...
%   'manyBlockers', 'funnel', 'fallingCells'
%'advectionInterpolation' like in interp2: 'spline', 'cubic', 'linear', 'nearest'
numCells = 99;
domainSize = 1; %m
simulationTime = 10; %s
fps = 120;

simulationSettings = struct(...
    'xCells', numCells, ...
    'yCells', numCells, ...
    'cellSize', domainSize / numCells, ...
    'numFrames', simulationTime * fps, ...
    'deltaTime', 1 / double(fps), ...
    'advectionAlgorithm', 'maccormack', ...
    'advectionInterpolation', 'linear', ...
    'dissipation', 1, ...
    'scenario', 'singleBlocker' ...
);

outputString = [simulationSettings.scenario, ' ', ...
    int2str(simulationSettings.xCells + 1), 'x', ...
    int2str(simulationSettings.yCells + 1), ...
    ' ', int2str(fps), 'fps'];

visualizationSettings = struct (...
    'fps', fps, ...
    'outputFolder', ['output\', outputString], ...
    'outputFile', [outputString '.avi'], ...
    'gridCellSize', 800 / numCells, ...
    'isobarWidth', 400 / numCells, ...
    'drawVelocity', false, ...
    'drawPressure', false, ...
    'distanceRefreshInterval', 1 ...
);

% clear old output
delete([visualizationSettings.outputFolder '\image*.png']);

disp(['Start simulation at: ', datestr(now)]);
[divergence, fluidPixels] = runSimulation(simulationSettings, visualizationSettings);
disp(['Finish simulation at: ', datestr(now)]);

% save and display profiling data
if enableProfiling
    profsave
end
% profile viewer

%% quality measures
figureFileName = [visualizationSettings.outputFolder '\' outputString];
figure('Name','Divergence Plot','NumberTitle','off')
plot(1:simulationSettings.numFrames - 1, divergence)
title('Sum of divergence')
xlabel('Frame')
ylabel('Divergence')
savefig([figureFileName ' divergence.fig'])
print([figureFileName ' divergence.png'],'-dpng')

figure('Name','Fluid Pixels Plot','NumberTitle','off')
plot(1:simulationSettings.numFrames - 1, fluidPixels ./ 1000)
title('Fluid Pixels *1000')
xlabel('Frame')
ylabel('Fluid Pixels')
savefig([visualizationSettings.outputFolder '\' outputString ' fluid pixels.fig'])
print([figureFileName ' fluid pixels.png'],'-dpng')
    
%% save resulting video
saveVideo(visualizationSettings, false);