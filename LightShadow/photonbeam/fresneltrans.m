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

function [T] = fresneltrans(cos_theta_i,eta2_over_eta1)
%Dielectric Fresnel transmission
%cos_theta_i: incom cos between normal and incoming direction
%eta2_over_eta1: relative index of refraction

    power = sign(cos_theta_i) * 2;
    c = abs(cos_theta_i);
    g_s = sqrt((eta2_over_eta1.^power)+c.*c-1);
    gmc = g_s-c;
    gpc = g_s+c;
    R = ((gmc./gpc).*(gmc./gpc))/2.*(1+(((c.*gpc-1)./(c.*gmc+1)).^2));
    R(imag(R)~=0)=1; 

    T = 1-R;
end

