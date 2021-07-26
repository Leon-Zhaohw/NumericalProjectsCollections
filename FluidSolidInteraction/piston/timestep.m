% ---------------------------
% filename timestep.m
%
% Used to compute the time step according to the
% Courant, Friedichs and Levy condition.
%
% it requires dxmin, the the minimum element size.
% wx is the vector of the velocity of the node
% vcelerity is the vector of sound celerity at each node
% vsol(:,2)./vsol(:,1) compute the vector of local fluid celerity
%                      at each node

function[Delta_t]=timestep(vcor,vsol,wx,vcelerity,CFL);

nnt     = size(vcor,1);
dxmin   = min(vcor([2:nnt],1)-vcor([1:nnt-1],1));
Delta_t = CFL.*dxmin./max(abs(vsol(:,2)./vsol(:,1))+vcelerity+abs(wx));
