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

function [ t, pdf_t ] = exponential_sampling(xi,sigma_t)
%xi: random number [0,1] 
%sigma_t: extinction coefficient

% Returns the position t of a sample along the extended source and its
% corresponding pdf according to an exponential pdf, pdf_t
    t = -log(1-xi)./sigma_t ; %equation (9) in [Habel13pbd]
    pdf_t = sigma_t  .* exp(-sigma_t.* t); %equation (9) in [Habel13pbd]
end