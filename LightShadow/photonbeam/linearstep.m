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

function [ out ] = linearstep( edge0,edge1,t )
%edge0: start of linear step function
%edge1: end of linear step function
%t: evalueation of ste function

    t = ((t - edge0)./(edge1 - edge0)); 
    t = max(min(t,1),0);

    out =  t;

end

