% PLOT THEORETICAL SINGLE SLIT PATTERN

function plotTheoreticalSingleSlit(N,a,b,twopi,x,x1)

% Plots the theoretical intenstity pattern for our single slit.
% Function is written by AIV. Version 15. October 2017

%figure;
theta = atan2(([1:N]-(N/2)),b);
betah = (twopi*a/2).*sin(theta);
sinbetah = sin(betah);
theoretical = (sinbetah./betah).*(sinbetah./betah);
x12 = x1(:,1).*x1(:,1);  % Calculate intensities
scaling = max(x12);
plot(x,theoretical.*scaling,'-g');
return;