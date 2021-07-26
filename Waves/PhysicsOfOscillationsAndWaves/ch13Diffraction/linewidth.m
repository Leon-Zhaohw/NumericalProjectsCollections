% CALCULATES LINEWIDTHS (FWHM) FOR SINGLE SLIT AND GAUSSIAN INTENSITY PROFILE

function  linewidth(N,lambda,x1)

% Calculates linewidths (FWHM) for single slit and Gaussian 
% intensity profile.
% Function is written by AIV. Version 15. October 2017

x12 = x1(:,1).*x1(:,1);  % Calculation of intensities

mx2 = max(x12(:,1))/2.0;
lower = 1;
upper = 1;
for k = 1:N-1
        if ((x12(k,1)<=mx2) && (x12(k+1,1)>=mx2)) 
            lower = k;
        end
        if ((x12(k,1)>=mx2) && (x12(k+1,1)<=mx2)) 
            upper = k;
        end
end
disp('FWHM: ')
(upper-lower)*1.0/lambda
return;