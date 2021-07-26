function [ frame, fluidPixels ] = calcFrame( visualizationSettings, ...
    distanceField, obstacle, velocity, pressureField )
%CALCFRAME Generate an image out of the distance field and obstacles using 
%       the visualization settings. If desired, the velocity and pressure 
%       field can be added as overlay.
%   return frame: matrix(yNodes * gridCellSize, xNodes * gridCellSize, colorChannel)
%          fluidPixels: scalar, number of fluid pixels
%   visualisationSettings: struct with these fields:
%       gridCellSize: inflate one cell to this many pixels
%       isobarWidth: width of the isobars (areas with similar distance)
%       drawVelocity: bool should display velocity field
%       drawPressure: bool should display pressure field
%   distanceField: matrix(yNodes, xNodes)
%   obstacle: boolean matrix(yNodes, xNodes)
%   velocity: matrix(yNodes, xNodes, dim)
%   pressureField: matrix(yNodes, xNodes)

    gridCellSize = visualizationSettings.gridCellSize;

    yNodes = size(distanceField, 1);
    xNodes = size(distanceField, 2);
    
    % fix errors in inputs
    distanceField(isnan(distanceField)) = 0;
    distanceField(isinf(distanceField)) = 0;
    pressureField(isnan(pressureField)) = 0;
    pressureField(isinf(pressureField)) = 0;
    
    [x, y] = meshgrid(1 : xNodes, 1 : yNodes);
    [xQuery, yQuery] = meshgrid(1 : 1 / gridCellSize : xNodes, ...
                               1 : 1 / gridCellSize : yNodes);
    inflatedDistanceField = interp2(x, y, distanceField, xQuery, yQuery);
    inflatedPressureField = interp2(x, y, pressureField, xQuery, yQuery);
    
    % get water filled areas
    waterBody = sign(inflatedDistanceField);
    waterBody(waterBody >= 0) = 1;
    waterBody(waterBody < 0) = 0;
    waterBody = logical(waterBody);
    
    xInflatedSize = size(inflatedDistanceField, 1);
    yInflatedSize = size(inflatedDistanceField, 2);
     
    % scale and round to get steps in 0...+inf
    distanceSteps = double(int32(abs(inflatedDistanceField))) ./ ...
        visualizationSettings.isobarWidth;
    % limit to n isobars
    maxDistance = 10 ./ visualizationSettings.isobarWidth;
    distanceSteps = max(min(distanceSteps, maxDistance), -maxDistance); 
    stepColorInfluenceFactor = 0.5;
    distanceSteps = distanceSteps * stepColorInfluenceFactor;

    waterColorFactors = [0.5 0.5 0.9];
    airColorFactors = [0.9 0.9 0.9];

    % calculate colors
    redFrame = distanceSteps;
    redFrame(waterBody == true) = waterColorFactors(1) ...
        .^ (redFrame(waterBody == true) + 1);
    redFrame(waterBody == false) = airColorFactors(1) ...
        .^ (redFrame(waterBody == false) + 1);

    greenFrame = distanceSteps;
    greenFrame(waterBody == true) = waterColorFactors(2) ...
        .^ (greenFrame(waterBody == true) + 1);
    greenFrame(waterBody == false) = airColorFactors(2) ...
        .^ (greenFrame(waterBody == false) + 1);

    blueFrame = distanceSteps;
    blueFrame(waterBody == true) = waterColorFactors(3) ...
        .^ (blueFrame(waterBody == true) + 1);
    blueFrame(waterBody == false) = airColorFactors(3) ...
        .^ (blueFrame(waterBody == false) + 1);

    % assemble color channels to image
    frame = zeros(xInflatedSize, yInflatedSize, 3);
    frame(:, :, 1) = redFrame;
    frame(:, :, 2) = greenFrame;
    frame(:, :, 3) = blueFrame;
    
    % draw pressure
    if visualizationSettings.drawPressure
        
        % normalize for each frame individually
        maxPressure = max(max(inflatedPressureField));
        minPressure = min(min(inflatedPressureField));
        
        frameAddition = inflatedPressureField;
        % to -1...+1
        frameAddition = (frameAddition - minPressure) / (maxPressure - minPressure);
        frameAddition = frameAddition * 2 + 1; % to 0...1
        frame(:, :, 1) = frame(:, :, 1) + frameAddition / 4;
    end
    
    % draw obstacle
    inflatedObstacle = double(obstacle);
    [xObst, yObst] = meshgrid(1 : xNodes, 1 : yNodes);
    [xObstQuery, yObstQuery] = meshgrid(1 : 1 / gridCellSize : xNodes, ...
                                        1 : 1 / gridCellSize : yNodes);
    inflatedObstacle = interp2(xObst, yObst, inflatedObstacle, ...
        xObstQuery, yObstQuery);
    inflatedObstacle = repmat(inflatedObstacle, [1 1 3]);
    frame(inflatedObstacle >= 0.5) = frame(inflatedObstacle >= 0.5) - 0.25;
    
    % draw velocity field
    % todo: find a way to not open a window for this
    % set(gcf,'Visible', 'off'); and figure('visible','off'); didn't work
    if visualizationSettings.drawVelocity
        frameXSize = size(inflatedDistanceField, 1);
        frameYSize = size(inflatedDistanceField, 2);
        drawAreaSize = 800; % should be displayable on every monitor
        fig = figure('Position',[0 0 drawAreaSize drawAreaSize], ...
            'Color', [0 0 0]);
        
        % set velocity at the borders to 0
        % should prevent the plot from being resized 
        % to fit arrows into the border
        velocity(1,:,:) = 0;
        velocity(yNodes,:,:) = 0;
        velocity(:,1,:) = 0;
        velocity(:,xNodes,:) = 0;
        
        quiver(x-1,y-1,velocity(:,:,1),velocity(:,:,2),0.5,'color',[1 1 1]);
        axis off
        set(gca,'position',[0 0 1 1],'units','normalized')
        F = getframe();
        velocityFieldFrame = double(F.cdata) / 255;
        velocityFieldFrame = imresize(velocityFieldFrame, ...
            [frameXSize frameYSize]);
        close
        
        velocityFieldFrame = flip(velocityFieldFrame, 1);

        if size(frame) == size(velocityFieldFrame) % window has minimum size
            frame = frame + velocityFieldFrame;
        end
    end
    
    fluidPixels = sum(sum(waterBody));
end

