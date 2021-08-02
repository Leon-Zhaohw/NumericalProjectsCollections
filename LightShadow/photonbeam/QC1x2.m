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

function res = QC1x2(eta)
%Polynomial approximations from [d'Eon and Irving 2011]
    if eta<1
        res = 0.919317 - 3.47930*eta + 6.753350*eta.^2 - 7.809890*eta.^3 ...
                                     + 4.985540*eta.^4 - 1.368810*eta.^5;
    else
        res = -9.23372 + 22.2272.*eta - 20.9292.*eta.^2 + 10.2291.*eta.^3 ...
                                     - 2.54396.*eta.^4 + 0.254913.*eta.^5 ;    
    end %equation (8) in [Habel13cid]
end