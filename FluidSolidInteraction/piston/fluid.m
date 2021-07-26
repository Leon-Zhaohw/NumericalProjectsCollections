% --------------------------------------------
% filename fluid.m
%
% we compute here the three fluid models.
% The pressure of the piston for each of them
% is stored in presL2t
%
% presL2t is a table with three columns.
% There are as many lines as simulation steps.
% Column 1 contains the pressure on the piston computed using fluid_model == 'A'
% Column 2 contains the pressure on the piston computed using fluid_model == 'B'
% Column 3 contains the pressure on the piston computed using fluid_model == 'C'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A-case: adiabaticity
% we compute and store in presL2t the pressure on the piston
% for the A model
%
presL2t(istep,1) = pres_init0*(L_0./(L_0+ u_t ))^(gam);

% We compute the total mass of the fluid for the A model
%
rho            = rho_init*(presL2t(istep,1)./pres_init)^(1./gam);
M_ACase(istep) = rho*A*(L_0+ u_t);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B-case
% we compute and store in presL2t the pressure on the piston
% for the B model
%
presL2t(istep,2)=pres_init(1)*(1+(gam-1)/2*(-u_dot_t./vcelerity(1)))^(2*gam/gamm1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C-case
%
% we update the mesh, using the mesh deformation technique

[wx,vcor_n,vcor_np1] = move_mesh(formulation,Delta_t,u_t-vsols0,u_dot_t,vsol,vcor_n,vcor_np1);
%
% If pure Eulerian on deformed mesh wx is set to 0.
% In this case wx is only used to move the mesh,
% it has no influence in the fluid model equations.
%
if (formulation==0)
    wx = wx.*0;
end
%
% Velocity of the piston from the fluid point of view
%
vsol(nnt,2) = vsol(nnt,1)*u_dot_t;
%
% Lax-Wendroff
% if ifct is set to zero we compute a solution of order 1
% if ifct is set to 2 we compute a solution of order 2
% if ifct is set to 1 we combine both low and high order solutions
%
ifct = 1;
[vsol, vflux] = Lax_Wendroff(gamm1,ifct,Delta_t,vsol,vpres,wx, conec, vcor_n, vcor_np1,number);
%
% We now update the vectors of temperature and pressure on each node of the mesh
%
vtemp            = temperature(C_v,vsol);
vpres            = pressure(R,vtemp,vsol);
%
% at this stage the corrected pressure on the piston is stored for the
% C model.
%
presL2t(istep,3) = vpres(nnt);
%
% We compute the total fluid energy variation (impulsion)
%
if(istep==1)
    Imp_fl(istep,:)=presL2t(istep,:)*A*u_dot_t*Delta_t;
else
    Imp_fl(istep,:)=Imp_fl(istep-1,:)+presL2t(istep,:)*A*u_dot_t*Delta_t;
end
%
% We store the data used for displaying and post-processing
% NB: all the data are not used for post-processing but they are
% stored if you want to use them.
%
histo_velocity(istep,[1 2])   = [u_dot_t wx(nnt)];
histo_pressure(istep,[1 2 3]) = presL2t(istep,:);
histo_deformation(istep,2)    = vcor_np1(nnt,1)-L_0;
%
% We compute now the global mass, global momentum, and global energy
%
M(istep)      = integ(vsol(:,1).*A,vcor_np1(:,1),conec);
QdM(istep)    = integ(vsol(:,2).*A,vcor_np1(:,1),conec);
NRJ_flu(istep)= integ(vsol(:,3).*A,vcor_n(:,1),conec);
