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

function res = QC2x3(eta)
%Polynomial approximations from [d'Eon and Irving 2011]
    if eta<1
        res = 0.828421 - 2.62051.*eta + 3.362310*eta.^2 - 1.952840.*eta.^3 ...
                                     + 0.236494.*eta.^4 + 0.145787.*eta.^5;    
    else
        res = -1641.1 + 135.926./eta.^3  - 656.175./eta.^2 + 1376.53./eta + 1213.67.*eta ... 
              -568.556.*eta.^2 + 164.798.*eta.^3 - 27.0181*eta.^4 + 1.91826*eta.^5;    
    end %equation (9) in [Habel13cid]
end