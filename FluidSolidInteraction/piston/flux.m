% -----------------------
% filename flux.m
%
% General formulation for a convective flux
%
% The section is directly included in the pressure term

function[foo]=flux(vsol,pres,wx);

foo = [vsol(:,2)-vsol(:,1).*wx  ...
    vsol(:,2).^2./vsol(:,1)+pres-vsol(:,2).*wx  ...
    (vsol(:,3)+pres).*vsol(:,2)./vsol(:,1)-vsol(:,3).*wx];

