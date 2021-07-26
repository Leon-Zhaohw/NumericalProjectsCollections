% ------------------------------
% filename: flu_mass.m
%
% used to compute the mass matrix

function [vme]=flu_mass(ie,vcore);

x   = vcore(2 ,:) - vcore(1,:);
xl  = sqrt(x * x' );

vme = xl./6.*[ 2 1
    1 2];
