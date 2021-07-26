% -----------------
% filename pressure.m
% used to compute the pressure at each node

function[vpres]=pressure(R,vtemp,vsol);

vpres = vsol(:,1).*R.*vtemp;
