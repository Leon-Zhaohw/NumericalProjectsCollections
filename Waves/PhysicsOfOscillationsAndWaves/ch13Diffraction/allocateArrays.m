% ALLOCATE ARRAYS WE NEED

function [x,x2,x0,x1,sc,r] = allocateArrays(nWavel,N)
% Allocates space for various arrays
% Function is written by AIV. Version 15. October 2017

x = linspace(-nWavel/2, nWavel/2, N); % A relative position 
							% array for plot
x2 = linspace(-N,N,2*N+1);  % Simil, but for plot/test of  
							% hjelp functions
x0 = zeros(N,2);            % Excitation data, amplitudes 
							% and phases
x1 = zeros(N,2);            % Amplitudes at screen, amplitudes 
							% and phases
sc = zeros(2*N + 1,2);    	% Store sin/cos for component 
							% calculations
r = zeros(2*N + 1,2);       % Distance-table: reduction 
							% factor and phase-correction
                            % based on path length
return;