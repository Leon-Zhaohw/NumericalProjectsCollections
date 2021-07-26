function antennaDiagram3

% Calculate the radiation pattern for a dipole antenna in log units
% in the radial direction.

N = 1024;
theta = linspace(-pi/2.0,3.0*pi/2.0,N);    % Angles
costheta = cos(theta);
intensities = costheta.*costheta;
intensities = log10(intensities*1000.0);
for i = 1:N
    if(intensities(i)<0)
        intensities(i)=0;
    end
end
polar(theta,intensities);