%******************************************************************************
%
%Copyright (c) 2013 by Disney-Pixar
%
%Permission is hereby granted to use this software solely for 
%non-commercial applications and purposes including academic or 
%industrial research, evaluation and not-for-profit media
%production.  All other rights are retained by Pixar.  For use 
%for or in connection with commercial applications and
%purposes, including without limitation in or in connection 
%with software products offered for sale or for-profit media
%production, please contact Pixar at tech-licensing@pixar.com.

%******************************************************************************

function [t, pdf] = equiangular_sampling(xi, u0, u1, h)
%xi:  random number [0,1]
%u0: start of ray
%u1: end of ray
%h: distance to ray of samplign point

% Equi-angular sampling along a ray with respect to a point as in e.g. [Kulla and Fajardo 2011].
% Returns the position t of a sample along the extended source and its
% corresponding pdf according to an equiangular pdf.

    % clamp h to avoid numerical singularity
    h = max(h, 1e-12);
	theta_min = atan(u0 ./ h);
	theta_max = atan(u1 ./ h);
    
	t = h .* tan(lerp(xi, theta_min, theta_max)); %equation (10) in [Habel13pbd]
	pdf = h ./ ((theta_max - theta_min) .* (h.^2 + t.^2)); %equation (10) in [Habel13pbd]

end