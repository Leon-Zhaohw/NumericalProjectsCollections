% CREATE WAVE PACKAGE IN SPACE (at t = 0)

function [x,z] = pg_wpack(N,xmax,xlambda,xsigma)

% Create a wave package in space (!). Version Oct 5 2017 AIV
% Input parameters: N: Points in the description, 
% xmax: Defines the interval x is defined ( |0,xmax>), 
% xlambda: spatial wavelength for the central wavelengh,  
% xsigma: the width in the gaussian shaped wave package.
% Returns: x array as well as the wave package array.

x = linspace(0,xmax*(N-1)/N,N);
xr = xmax/8.0; % Startpoint for the center of the wave package
xfreq = 1/xlambda;  % Spatial frequency
y = cos((x-xr)*2*pi*xfreq);
convol = exp(-((x-xr)/xsigma).*((x-xr)/xsigma));
z = y.*convol;
return;
