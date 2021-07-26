% ------------------------------
% filename spring.m
%
% we compute displacement u_t, velocity u_dot_t and acceleration u_double_dot_t
% of the piston
%
function[u_t, u_dot_t, u_double_dot_t]=spring(vprel,Ppiston,pres_init0,A,Delta_t,u_t,u_dot_t,u_double_dot_t);

%============================================================
% incremental scheme used:
% KT    = 4M +Delta_t^2K,
% KT.DU = Delta_t^2.F+4MDelta_t.VIT+Delta_t^2M.ACC-Delta_t^2K.VSOL
%============================================================

vkg  = vprel(1); % spring rigidity
vmg  = vprel(2); % mass of the piston
vfg  = (Ppiston-pres_init0).*A;

% -----
% two intermediate values for computing purpose
vkt  = 4.*vmg + Delta_t^2*vkg;
vres = Delta_t^2*vfg'+4*Delta_t*vmg*u_dot_t+Delta_t^2*vmg*u_double_dot_t-Delta_t^2*vkg*u_t;

%----- resolution
% We compute the variation of position between time step n and n+1
vdu = vkt\vres;
%
% we update u_t, u_double_dot_t and u_dot_t
%
u_t            = u_t + vdu;        %
inter          = 4./Delta_t^2*vdu-4./Delta_t*u_dot_t-u_double_dot_t;
u_dot_t        = u_dot_t+Delta_t./2*(u_double_dot_t+inter);
u_double_dot_t = inter;



