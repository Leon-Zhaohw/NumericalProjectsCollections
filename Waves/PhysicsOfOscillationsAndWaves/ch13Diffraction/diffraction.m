function diffraction(code)  % MAIN PROGRAM

% This program calculates and plots intensity patterns for a 
% variety of diffraction and/or interpherence phenomena with 
% cylindrical symmetry.
% Functionalities: code = 1: One slit, 2: Gaussian intensity 
% profile, 3: Straight edge, 4: Double slit, 5: Read excitation 
% data from file (amplitude + phase)
% Program is written by AIV. Version 15. October 2017

% Establishes essential parameters for the calculations.
% Results depend critically on these. See the code for this 
% function. 

[lambda,a,b,nWavel,N,twopi,Nhalf] = parameters;

% Allocates arrays for the calculations:
[x,x2,x0,x1,sc,r] = allocateArrays(nWavel,N);

% Generate or read in excitation data:
[x0] = generateExcitation(code,lambda,a,N,Nhalf,twopi,x0);

% Calculates sines, cosines, distances and relative phase 
% differences for vectors between the plane of excitation and 
% the screen for the final pattern:
[sc,r] = generateRelPositionData(N,b,lambda,twopi);

% Sum all contributions to every point on the screen (main 
% loop):
[x1] = summation(N,x0,r);

% Plots intensities for diffraction pattern along with a 
% marking of the excitation
plotDiffraction(x,x0,x1);

% Calculates and write out linewidths in case the excitation 
% was a single slit or Gaussian profile, and write to 
% screen the actual linewidth of the intensity profile.
if (code==1) || (code==2) 
    linewidth(N,lambda,x1);
end

% Plots expected theoretical intensity profile for a single 
% slit:
if code==1
    plotTheoreticalSingleSlit(N,a,b,twopi,x,x1);
end

% Option: Save data to a file (as a string of floating point 
% numbers):
%writeToFile(x1);

% Removes all plots when we leave the program (cleans up):
input('Close all figures');
close all