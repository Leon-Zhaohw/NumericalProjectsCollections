% ------------------------------
% filename: Lax_Wendroff.m
%
% We recall that vsol:
% column 1: the volumic mass
% column 2: the momentum
% column 3: the total energy (per volum unit)

function [vsol, vflux]=Lax_Wendroff(gamm1,ifct,Delta_t,vsol,vpres,wx, conec, vcor_n, vcor_np1,number);

[nnt, ndln]  = size(vsol);
[nelt, nnel] = size(conec);
%
vres = zeros(nnt,ndln);
%
% we compute the nodal fluxes
%
% MANU voir flux.m
%
vflux  = flux(vsol,vpres,1.*wx);
%
% Loop over the elements: 2 steps are included
%
vmg_n=zeros(nnt);
vmg_np1=zeros(nnt);
%
for ie=1:nelt

    kloce      = conec(ie,1:nnel);
    vcore_n    = vcor_n(kloce,1);
    vcore_np1  = vcor_np1(kloce,1);
    vflue      = vflux(kloce,:);
    wxe        = wx(kloce);
    %
    % we compute the elementary residual
    %
    vrese         = flu_residual(ie,Delta_t,conec,vsol,vpres,vcore_n,vflue,wxe,gamm1);
    %
    % Then we proceed to the assembling
    %
    vres(kloce,:) = vres(kloce,:)+vrese;
    %
    % Mass matrix at time station 'n'
    %
    vmg_n(kloce,kloce) = vmg_n(kloce,kloce)+flu_mass(ie,vcore_n);
    %
    % Mass matrix at time station 'n+1'
    %
    vmg_np1(kloce,kloce) = vmg_np1(kloce,kloce)+flu_mass(ie,vcore_np1);
    %
    % We have to take care of the boundary conditions for the fluxes terms
    % for the first finite element (ie=1) and the last one (ie = nelt)
    %
    if ie==1    vres(1,:)  = vres(1,:)  + Delta_t.*vflue(1,:);  end
    if ie==nelt vres(nnt,:)= vres(nnt,:)- Delta_t.*vflue(2,:);  end
end		% end of the loop over the elements

%
% 'diagonalization' of the mass matrix
% for the iterative resolution scheme
%
xlumpm = sum(vmg_n)';
%
% We compute the increment of the solution with a shock capturing technique
%
du = shock_capture(ifct,vcor_n,conec,vmg_n,vmg_np1,vres,vsol,xlumpm,number);
%
% We take into account the boundary condition
%
du(1,2)   = 0.; % The velocity on node number 1 (closed wall) set to 0.0
du(nnt,2) = 0.; % The velocity on node number nnt has been already updated by the struture solver.
%
% At this stage we update the solution (volumic mass, momentum, total energy)
%
vsol = vsol + du;
