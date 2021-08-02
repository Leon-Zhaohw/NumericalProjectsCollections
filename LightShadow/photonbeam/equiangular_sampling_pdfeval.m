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

function pdf = equiangular_sampling_pdfeval(t, u0, u1, h)
%t: distance along extended source   
%u0: start of ray
%u1: end of ray 
%h: distance to ray of samplign point
% Evaluates the equiangular pdf for a given distance t along the extended source .   

    h = max(h, 1e-12);
	theta_min = atan(u0 ./ h);
	theta_max = atan(u1 ./ h);
	pdf = h ./ ((theta_max - theta_min) .* (h.^2 + t.^2)); %equation (10) in [Habel13pbd]
    
end