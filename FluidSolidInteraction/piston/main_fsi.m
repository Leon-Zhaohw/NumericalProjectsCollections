%-----------------------------------------------------------------------------
% This file main_fsi.m is the main script provided with the paper entitled
% 'An Initiation to fluid-structure interaction: application to the piston problem'
%
% 'fsi' is an abbreviation for 'fluid structure interaction'.
%
% It is a script coupling a structure solver and a fluid solver
% chosen in a list of three fluid solvers.
%
% The textbook we consider:
% a gas contained in a one dimension chamber, closed on its
% right part by a moving piston and closed on its left by a fixed wall.
% The piston is attached to an external and fixed point with a spring.
%
%            the one dimension chamber
%        |---------------------------------
%        |                                 _
% fixed  |                                | | the piston
% wall   |                                | |
%        |      enclosed air              | |----/\/\/\/\/\/\--o the fixed point
%        |                                | |    the spring
%        |                                | |
%        |                                 -
%        |---------------------------------
%
% The three fluid models:
%   o A-model perfect gas law with adiabaticity condition: P.V^gamma = cst
%   o B-model piston analogy model
%   o C-model 1D compressible fluid flow evolution with three fluid flow formulations:
%               0) Eulerian
%               1) Lagrangian
%               2) ALE (Arbitrary Lagrangian Eulerian)
%     Note that a shock capturing technique is used to ensure
%     fluid spatial stability for the C-model.
%
% -------------------------------
% clean the environment first
clear all;
close all;
% -------------------------------
% Initialization of the raw data:
%
% fluid first, structure second.
initialize_data_fluid;
initialize_data_structure;
%
Tmax  = 0.5*T0;  % We set the maximum simulation time.
%                  T0, the natural period of the piston
%                  is set in "initialize_data_structure.m" script.
CFL   = 0.8;	 % We set the CFL criterion (Courant, Friedichs and Levy)
%
% note: to make simple, at each time period all the fluid models are computed.
%       See 'fluid.m' script. The character variable 'fluid_model' permits to chose
%       the computed pressure used for the structure computation.
%       See 'structure.m' script.
%
% We choose here the fluid model.
% Comment and uncomment lines as you whish.
%
%fluid_model = 'A'; % A-model chosen.
%fluid_model = 'B'; % B-model chosen.
fluid_model = 'C'; % C-model chosen.
%
% -----------------
% We choose here the formulation when fluid_model = 'C'.
% Comment and uncomment lines as you whish.
%
% formulation = 0 % Eulerian with mesh deformation technique
% formulation = 1 % Lagrangian
formulation = 2;  % ALE (Arbitrary Lagrangian Eulerian)
%
% Note: to be coherent, formulation>0 involves fluid_model = 'C'.
%       The mesh deformation technique is computed whatever is the model
%       but this has no impact on the computation of the pressure for
%       fluid_mode = 'A' or fluid_model = 'B'.
%       Indeed, if (formulation=0) then fluid_model = 'A' or fluid_model =
%       'B', so p(t) the pressure on the piston only depends on its position or
%       velocity.

if (formulation>0)
    fluid_model = 'C';
end

i3Dview=1;
if(i3Dview==0) fprintf(1,'3D view INACTIVE! \n'); end
%
fprintf(1,'Fluid model = %c  \n',fluid_model);
fprintf(1,'Fluid formulation (0=Eulerian, 1=Lagrangian, 2=ALE) = %g  \n',formulation);

% -------------------
%
% We compute the initial time step, Delta_t, using the CFL criterion
% we also compute the minimal segment size dxmin.
% Note 1: initialy all the discretization points have a celerity wx=0,
% Note 2: since the mesh changes, Delta_t has to be computed at each iteration
%         of the simulation loop (see below).
% vcor_np1 coordinates of discretization points at time station 'n+1'
% vcor_np  coordinates of discretization points at time station 'n'
% initially they are equal
dxmin   = min(vcor_np1([2:nnt],1)-vcor_np1([1:nnt-1],1));
Delta_t = CFL.*dxmin./max(abs(vsol(:,2)./vsol(:,1))+vcelerity);

%
% We compute this modular for managing 3D display
% We display one result over ten.
modulo = ceil(ceil(Tmax./Delta_t)./10);
%
% Initialization of indices and time to zero
istep      = 0;  % iteration step counter
Total_time = 0.; % total simulation time
its        = 0;  % just for managing storage.
% --------------------------------
% begin of the simulation loop
while (Total_time<( Tmax - Delta_t))
    istep = istep + 1;
    %
    % we compute Delta_t according to CFL conditions
    Delta_t                = timestep(vcor_np1,vsol,wx,vcelerity,CFL);
    Delta_t_storage(istep) = Delta_t;
    Total_time             = Total_time + Delta_t;
    t(istep)               = Total_time;
    %
    %%%%%%%%%%%% STRUCTURE %%%%%%%%%%%%
    structure;
    %
    %%%%%%%%%%%% FLUID %%%%%%%%%%%%
    fluid;
    %
    %%%%%%%%%%%% Displaying results %%%%%%%%%%%%
    %
    fsi_display;
end
% end of the simulation loop
% ---------------------------------
% post-treatment, displaying results.
post_treatment;
