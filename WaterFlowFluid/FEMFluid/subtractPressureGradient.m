function [ newVelocityField ] = subtractPressureGradient( ...
    velocityField, pressure, obstacle )
%SUBTRACTPRESSUREGRADIENT Apply pressure on the velocity field
% based on: https://github.com/candycat1992/2DFluidSim
%   return newVelocityField: matrix(xNodes, yNodes, dim)
%   velocityField: matrix(yNodes, xNodes, dim)
%   pressure: matrix(yNodes, xNodes)
%   obstacle: boolean matrix(yNodes, xNodes)

    gradient = calcGradient( pressure, obstacle );
    
    % subract gradient from velocity for new velocity
    newVelocityField = velocityField - gradient;
	
end