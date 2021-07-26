% CREATE A SUM OF HARMONIC SPATIAL WAVES AT A GIVEN TIME t}

function [zrecon] = pg_wave(x,t,N,A,phase,k,omega,imin,imax)

% Generate the complete spatial wave using the Fourier 
% coefficients. Version Oct 5 2017 AIV
% Input parameters: x: position array, t: current time, 
% N: number of points, [A, phase, k]: amplitude, phase and 
% wavenumber arrays, respectively, omega: the dispersion 
% relation omega(k), [imin, imax]: minimum andmaximum index 
% that will be used in the arrays A, phase and k.
% Returns: zrecon: the position of the marker which give the
% position to where a peak with the central wavelength would 
% have ended up (for verification of proper functioning).
    
zrecon = zeros(1,N);
for i = imin:imax  % Sum over Fourier elements
    arg = k(i)*x - omega(i)*t + phase(i);
    zrecon = zrecon + A(i)*cos(arg);
end
return;