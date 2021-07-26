% FREQUENCY ANALYSIS OF WAVE PACKAGE IN SPACE}

function [A,theta,k] = pg_fft(z,N,xmax)

% Frequency analysis of a wave package in space. 
% Version Oct 5 2017 AIV
% Input parameters: z: the array describing the wave package,
% N: number point in this description, xmax: describe the x 
% interval |0, xmax>. Returns: Amplitude A and phase (theta) 
% in the frequency analysis as a functon of the wavenumber k.

Zf = fft(z)/N;
A = 2.0*abs(Zf);  % Ignore the error in Zf(1), don't use it
theta = atan2(imag(Zf),real(Zf));
xsamplf = N/xmax;   % Spatial sampling frequency
xfreq = linspace(0,xsamplf*(N-1)/N,N); % Spatial frequency
k = zeros(1,N);
k = 2.0*pi*xfreq;

% NOTE: Use the reminder of this function when you need to 
% pick frequency components for your wave package. You need 
% this in order to choose imin and imax in the program pg3.m
%figure;
%plot(A,'.-r'); % Plot to be able to choose points to be used
%plot(xfreq,A,'.-r');  % Alternative plot
%plot(xfreq,fase,'.-k');
return;