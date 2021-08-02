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


%This code is an example implementation of "Photon Beam Diffusion: A Hybrid
%Monte Carlo Method for Subsurface Scattering" [Habel13pbd]
%http://zurich.disneyresearch.com/~wjarosz/publications/habel13pbd.html

%For a description of the searchlight problem and the theory of diffusion 
%in context of light scattering, please consult:
%"Classical and Improved Diffusion Theory for Subsurface Scattering" [Habel13cid]
%http://zurich.disneyresearch.com/~wjarosz/publications/habel13cid.html

%References to equations in both publications are given throughout the code.

%The example code generates BSSRDF profiles for the searchlight problem and
%compares them to the Monte Carlo reference solution.
%The Monte Carlo solution has been generated according to MCML 
%http://omlc.ogi.edu/software/mc/
%http://www.atomic.physics.lu.se/biophotonics/our_research/monte_carlo_simulations/gpu_monte_carlo/
 
%Set scattering parameters. This is the set for the included Monte Carlo profile
sigma_a = 0.1;
sigma_s = 1;
g = 0;
eta = 1.33;

%Generate an dense array of radial distances to evaluate the searchlight problem
%for different BSSRDFs
r_array = 0.002:0.004:8;
r_array = r_array';

%Evaluate diffuse reflectance of the single scattering
%Single scattering is exact
single_p = single_scattering_profile(sigma_a,sigma_s,g,eta,r_array);

%To compare to the Monte Carlo solution, we need to take the Fresnel
%boundary into account for the incident orthogonal light. So we multiply by a
%factor of
Li = fresneltrans(cos(0),eta/1);


%Evaluate classical dipole profile
classical_p = Li.*classical_dipole(sigma_a,sigma_s,g,eta,r_array);

%Evaluate improved dipole profile
improved_p = Li.*improved_dipole(sigma_a,sigma_s,g,eta,r_array);

%Evaluate photon beam diffusion profile
PBD_p = Li.*PBD_profile(sigma_a,sigma_s,g,eta,r_array);


%Contains Monte Carlo reference for the given scattering parameters
load('MonteCarlo_0.1.mat');

%Plot multiple scattering and multi+single scattering
%in comparison to the Monte Carlo reference.
%Note that the plot is over a larger range than in [Habel13pbd]
fig1 = figure(1);
set(fig1, 'Position', [50 50 1200 800], 'Name', 'BSSRDF comparison plots');

subplot(2,2,1);
semilogy(r_array,r_array.*MC_p,'Color','black');
hold all;
semilogy(r_array,r_array.*classical_p,'Color','red');
semilogy(r_array,r_array.*improved_p,'Color',[0.2 0.2 0.8]);
semilogy(r_array,r_array.*PBD_p,'Color',[0.2 0.8 0.2]);

axis([0 8 10^-4 5*10^-2])
xlabel('x');
ylabel('x * R(x)');
legend('Monte Carlo','Classical Dipole','Improved Dipole','Photon Beam Diffusion');
title(['Multi Scattering: \sigma_s = ' num2str(sigma_a) ' \sigma_a = ' num2str(sigma_a) ', g=' num2str(g)], 'FontWeight', 'bold');

subplot(2,2,2);
semilogy(r_array,r_array.*(MC_p+single_p),'Color','black');
hold all;
semilogy(r_array,r_array.*(classical_p+single_p),'Color','red');
semilogy(r_array,r_array.*(improved_p+single_p),'Color',[0.2 0.2 0.8]);
semilogy(r_array,r_array.*(PBD_p+single_p),'Color',[0.2 0.8 0.2]);

axis([0 8 10^-4 10^-1])
xlabel('x');
ylabel('x * R(x)');
legend('Monte Carlo','Classical Dipole','Improved Dipole','Photon Beam Diffusion');
title(['Multi+Single Scattering: \sigma_s = '  num2str(sigma_a) ' \sigma_a = ' num2str(sigma_a) ', g=' num2str(g)], 'FontWeight', 'bold');

subplot(2,2,3);
semilogy(r_array,r_array.*MC_p,'Color','black');
hold all;
semilogy(r_array,r_array.*classical_p,'Color','red');
semilogy(r_array,r_array.*improved_p,'Color',[0.2 0.2 0.8]);
semilogy(r_array,r_array.*PBD_p,'Color',[0.2 0.8 0.2]);

axis([0 3 5*10^-3 5*10^-2])
xlabel('x');
ylabel('x * R(x)');
legend('Monte Carlo','Classical Dipole','Improved Dipole','Photon Beam Diffusion');
title(['Multi Scattering: \sigma_s = ' num2str(sigma_a) ' \sigma_a = ' num2str(sigma_a) ', g=' num2str(g) ' closeup'], 'FontWeight', 'bold');

subplot(2,2,4);
semilogy(r_array,r_array.*(MC_p+single_p),'Color','black');
hold all;
semilogy(r_array,r_array.*(classical_p+single_p),'Color','red');
semilogy(r_array,r_array.*(improved_p+single_p),'Color',[0.2 0.2 0.8]);
semilogy(r_array,r_array.*(PBD_p+single_p),'Color',[0.2 0.8 0.2]);

axis([0 3 5*10^-3 10^-1])
xlabel('x');
ylabel('x * R(x)');
legend('Monte Carlo','Classical Dipole','Improved Dipole','Photon Beam Diffusion');
title(['Multi+Single Scattering: \sigma_s = '  num2str(sigma_a) ' \sigma_a = ' num2str(sigma_a) ', g=' num2str(g) ' closeup'], 'FontWeight', 'bold');




