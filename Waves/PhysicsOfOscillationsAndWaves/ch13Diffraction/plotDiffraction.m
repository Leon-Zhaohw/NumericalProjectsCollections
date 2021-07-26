% PLOT THE DIFFRACTION PATTERN

function plotDiffraction(x,x0,x1)

% Plots intensities for diffraction picture along with a 
% marking of the excitation. Some extra possibilities are 
% given, for testing or special purposes.
% Function is written by AIV. Version 15. October 2017

%plot(x,x1(:,1),'-r');   % Plots amplitudes (red) (can 
						 % often be skipped)
figure;
x12 = x1(:,1).*x1(:,1);  % Calculation of intensities
hold on;
scaling = (max(x12)/8.0);
plot(x,x0(:,1).*scaling,'-r');   % Plot initial excitaion 
plot(x,x12(:,1),'-b');  % Plot relative intensities (blue)
xlabel('Position on screen (given as # wavelengths)');
ylabel('Relative intensities in the diffraction pattern');

% figure;
% plot(x,x1(:,2),'-k');   % Plot phases (black) (most 
						  % often skipped)
return;	