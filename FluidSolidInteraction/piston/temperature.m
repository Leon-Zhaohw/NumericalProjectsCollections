% ------------------------
% filename: temperature.m
%
%  We compute the temperature at each node
% vsol :
% column 1: the volumic mass
% column 2: the momentum
% column 3: the total energy (per volum unit)
% C_v     :the specific heat capacity of the gas

function [vtemp]=temperature(C_v,vsol);

vtemp=(vsol(:,3)./vsol(:,1)-0.5.*vsol(:,2).^2./vsol(:,1).^2)./C_v;
