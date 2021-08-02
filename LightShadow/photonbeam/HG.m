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

function [ p_HG ] = HG( g, costheta )
%g: mean cosine
%costheta: cos to direction of the ray

%Henyey-Greenstein phase function

    p_HG = 1./(4.*pi).* (1-g.*g)./((1+g.*g-2.*g.*costheta).^(1.5));

end

