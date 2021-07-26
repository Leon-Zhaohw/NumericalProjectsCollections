% CHOOSE PARAMETERS FOR THE CALCULATIONS

function [lambda,a,b,nWavel,N,twopi,Nhalf] = parameters
 
% Choose parameters for the calculation. 
% Written by AIV. Version 15. October 2017

% Choose resolution, distance to screen, and the width of the 
% area on screen the calculation should include. Some constants 
% are defined. Results depend critically on the parameters 
% set by this function. 
% Whether the result will be mainly a Fresnel- og Franuhofer 
% diffraction depend on the b parameter. nWavel must be 
% increased if b is large to include the full diffraction 
% pattern within the calculated area on screen.
% The parameters given in this particular code is suitable 
% for a double slit in the Fresnel regime (quite complicated 
% pattern). 

lambda = 4;         % Four points per wavelength resolution in 
					% excitation points
a = 20;             % Width of single slit, given in 
					% # wavelengths
b = 4000 * lambda;  % Distance to screen is b wavelengths
nWavel = 1024*3/2;  % # wavelengths along the screen (an 
					% integer!)
N = nWavel*lambda;  % Width of excitation area as well as 
					% screen in # wavelengths 
twopi = 2.0*pi;     % Somewhat innecessary, but speeds up 
					% a bit...
Nhalf = N/2;
return;