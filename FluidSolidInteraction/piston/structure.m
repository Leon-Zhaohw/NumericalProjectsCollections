% --------------------------------------
% filename structure.m
%
%
vsols0 = u_t;
%
% note: to make simple, at each time period all the fluid models are computed.
%       We state a correspondance between the 'fluid_model'
%       and 'indice_fluid_model' variables.
%
% The table Ppiston contains the three computed pressure at each time station.
if (fluid_model == 'A')
    indice_fluid_model = 1;
end
if (fluid_model == 'B')
    indice_fluid_model = 2;
end
if (fluid_model == 'C')
    indice_fluid_model = 3;
end
%
% presL2t is a table with three columns.
% There is as many line as simulation steps.
% Column 1 contains the pressure on the piston computed using fluid_model == 'A'
% Column 2 contains the pressure on the piston computed using fluid_model == 'B'
% Column 2 contains the pressure on the piston computed using fluid_model == 'C'
% presL2t is updated in the fluid.m script.

if(istep==1) % particular case when we start the simulation
    Ppiston = presL2t(1,indice_fluid_model);
else
    Ppiston = presL2t(istep-1,indice_fluid_model);
end
%
% we compute:
% displacement u_t, velocity u_dot_t and acceleration u_double_dot_t
% of the piston
[u_t, u_dot_t, u_double_dot_t]=spring(vprel,Ppiston,pres_init0,A,Delta_t,vsols0,u_dot_t,u_double_dot_t);
%
% We store the displacement for plotting purpose
%
histo_deformation(istep,1) = u_t;
%
% We compute the mechanical energy of the structure
% (kinetic part and potential part)
%
Ec(istep) = 1/2*vprel(2)*u_dot_t^2;
Ep(istep) = 1/2*vprel(1)*(Lspe-u_t-Lsp0)^2;  
Em(istep) = Ec(istep)+Ep(istep);
%
% We compute the external force and the power exerted on the piston
% =-[p(x)*A]_0^L(t)
% =-[p(x)*A*v(x)]_0^L(t)
Force_ext(istep)     = -(vpres(nnt)-vpres(1)).*A;
