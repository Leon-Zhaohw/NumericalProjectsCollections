% ----------------------
% filename flu_residual.m
% 
function [vrese]=flu_residual(ie,Delta_t,conec,vsol,vpres,vcore0,vflue,wxe,gamm1);

% Extraction of the parameters
%
kloce  = conec(ie,[1 2]); % connectivity of the element
vsole  = vsol(kloce,:);   % solution on the element
vprese = vpres(kloce,:);  % pressure on the element

x = vcore0(2 ,:)-vcore0(1,:); xl0=sqrt(x * x' );
%
% step 1 elementary: we compute the intermediate solution
% at t + Delta_t/2 on the element
%
vsol1_2=     0.5.*(vsole(1,:)+vsole(2,:)) ...
                -Delta_t.*0.5./xl0.*(vflue(2,:)-vflue(1,:)) ;
%
% We compute the nodal value using the average values on the elements
%
vsolmoy = (vsole(1,:)+vsole(2,:)).*0.5;
vsnod   = vsole + [vsol1_2-vsolmoy ; vsol1_2-vsolmoy];
%
% we compute the pressure (the section A is included)
%
vpres1_2S = gamm1.*(vsnod([1:2],3)- ...
                    vsnod([1:2],2).^2.*0.5./vsnod([1:2],1));

%
% We compute the intermediate flux
%
vflu1_2 = flux(vsnod,vpres1_2S,wxe);
%
% step 2 (global): we compute the elementary residual
% 
vrese = Delta_t.*0.5.*[-1 -1 ; 1 1 ]*[vflu1_2(1,:); vflu1_2(2,:)];

