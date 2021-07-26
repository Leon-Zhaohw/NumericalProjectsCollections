% SUMMATION OF ALL CONTRIBUTIONS

function [x1] = summation(N,x0,r)
% Runs through x1 (screen) from start to end and sum all 
% contributions from x0 (the excitation line) with proper 
% amplitude and phase.
% Function is written by AIV. Version 15. October 2017
    
for n = 1:N    
    relPos1 = N+2-n;
    relPos2 = relPos1+N-1;
    amplitude = x0(:,1).*r(relPos1:relPos2,1);
    fase = x0(:,2) - r(relPos1:relPos2,2);
    fasor(:,1) = amplitude .* cos(fase);
    fasor(:,2) = amplitude .* sin(fase);
    fasorx = sum(fasor(:,1));
    fasory = sum(fasor(:,2));
    x1(n,1) = sqrt(fasorx*fasorx + fasory*fasory);
    x1(n,2) = atan2(fasory, fasorx);   
end
return;