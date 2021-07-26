% ============================================================
% filename: initialize_data_fluid.m
%
% We define here: the geometry, the mesh, the fluid properties,
%                 the initial conditions, etc.
%
% We call here temperature.m and pressure.m scripts
%===========================================================
%
% we set the data size for the C-model
nelt = 70 ;       % nelt number of finite elements for C-model
nnt  = nelt+1;     % nnt  number of nodes (i.e. discretization points)
nnel = 2;          % nnel number of node per finite element
ndln = 3;          % ndln number of dof per node
ndle = ndln.*nnel; % ndle number of dof per element
ndlt = nnt.*ndln;  % ndlt total number of dof
%
% -------------------------
% Geometry of the chamber
%
L_0 = 1;  % length of the chamber at rest (meters)
A   = 1;  % piston section

% --------------------------
% Note: the initial condition may be changed here
U0   = 0.2;       % initial condition, i.e. initial displacement of the piston
L_t  = L_0 + U0;  % length of the chamber under initial condition

% -------------------------
%
% we set the initial coordinates of the discretization points
x        = [0:L_t/(nnt-1):L_t];
%
vcor     = x';    % vector of coordinates
vcor0    = vcor;  % saving this value for plotting purposes (see fsi_display.m)
vcor_n   = vcor;  % coordinates of discretization points at time station 'n'
vcor_np1 = vcor;  % coordinates of discretization points at time station 'n+1'
% the connectivity table :
% e.g: conec(1,:) = [1,2], connec(2,:) = [2,3], etc.
%      meaning that node number 1 and 2 belong to finite element number 1
%      and node number 2 and 3 belong to finite element number 2, etc.
conec([1:1:nelt],1)=[1:1:(nnt-1)]';
conec([1:1:nelt],2)=[2:1:nnt]';

%
% number(i)=k, where i is a node number,
% k is the number of finite element it belongs to.
%
number = zeros(nnt,1);
for i=1:nnt
    number(i) = length(find(conec==i));
end
%
% QUESTION MANU quid ?
iopen1 = 0; % Node number 1 is on the left of the chamber.
% It is a priori connected to the fixed wall.
% Except for fluid_model = 'B'.
% The variable iopen1 is used to distinguish this case.
% iopen1 = 0 corresponds to fluid_model = 'A' and fluid_model = 'C'
% iopen1 = 1 corresponds to fluid_model = 'B'.
%============================================================
% Physical properties of the flow
%===========================================================
%
gam   = 1.4;     % the specific heat ratio of the gas
gamm1 = gam-1.;
R     = 287;     % the individual gas constant
C_v   = R/gamm1; % the specific heat capacity of the gas

pres_init0 = 1E5;  % initial pressure for chamber length = L0
temp_init0 = 300;  % initial temperature
rho_init0  = pres_init0/gamm1/C_v/temp_init0; % initial volumic mass
%
% we compute the initial pressure at rest using PV^gamma=cste
% we compute volumic mass and temperature at rest.
% We set the external pressure.
if(iopen1 == 0)
    pres_init = pres_init0*(L_0./L_t).^gam;
end

rho_init  = rho_init0*(pres_init/pres_init0).^(1./gam);
temp_init = pres_init./rho_init/gamm1/C_v;
p_ext     = 0*pres_init0; % pressure on the right of the piston
%
% we set the initial fluid velocity
% and the initial total fluid energy
u_init = 0.;
e_init = pres_init./gamm1./rho_init+0.5*u_init.^2.;

% --------------------------------------------------
% Initialization of the table vsol used to store the solutions
% on each node.
% column 1: the volumic mass
% column 2: the momentum
% column 3: the total energy (per volum unit)
%
% Note that there are as many columns for vsol as ndln the number of dof per node
%
vsol      = zeros(nnt,3);
vsol(:,1) = rho_init.*ones(nnt,1);
vsol(:,2) = rho_init.*u_init.*ones(nnt,1);
vsol(:,3) = rho_init.*e_init.*ones(nnt,1);

% ------------------------------
% We set vectors for temperature, pressure and sound celerity
% at each node.
vtemp     = temperature(C_v,vsol);
vpres     = pressure(R,vtemp,vsol);
vcelerity = sqrt(gam.*gamm1.*C_v.*vtemp);
cel_init  = vcelerity(1);

%============================================================
% Initialization
%===========================================================
%
% vflux: convective fluxes vector
% wx: vector of the node velocity
% presL2t: table for storing the piston pressure computed by each fluid model
%
vflux               = zeros(nnt,ndln);
wx                  = zeros(nnt,1);
presL2t(1,[1 2 3])  = pres_init.*[1 1 1];
