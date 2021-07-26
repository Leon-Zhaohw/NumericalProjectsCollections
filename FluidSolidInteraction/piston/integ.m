% ----------------
% filename integ.m
%
% Used to compute the numerical integration.
% We perform the numerical integration according to the finite element method.
%
% This permits to compute the global mass, global momentum, and global energy
%
%
function[foo]=integ(vsol,x,conec);

% nelt number of finite elements
%
nelt = size(conec,1);
foo  = 0.;
%
% Loop over the elements.
%
for ie=1:nelt
    kloce = conec(ie,:);
    vcore = x(kloce);
    xl    = abs(vcore(2)-vcore(1));
    vsole = vsol(kloce);
    foo   = foo+vsole'*[1;1]*xl./2;
end
