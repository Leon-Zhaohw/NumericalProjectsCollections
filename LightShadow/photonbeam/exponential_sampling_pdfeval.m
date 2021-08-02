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

function [ pdf ] = exponential_sampling_pdfeval(t,sigma_t)
%t: distance along extended source
%sigma_t: extinction coefficient

% Evaluates the exponential pdf for a given distance along the extended source    
    pdf =  sigma_t .* exp(- sigma_t.* t);   %equation (9) in [Habel13pbd]
end