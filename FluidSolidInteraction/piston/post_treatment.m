% ---------------------------------------
% filename post_treatment.m
%
% note: to make simple, at each time period all the fluid models are computed.
%       We state a correspondance between the 'fluid_model'
%       and 'indice_fluid_model' variables.
%
if (fluid_model == 'A')
    indice_fluid_model = 1;
end
if (fluid_model == 'B')
    indice_fluid_model = 2;
end
if (fluid_model == 'C')
    indice_fluid_model = 3;
end

npas = istep;
%
% normalisation of force, energy and mass
%
F0      = vprel(1)*(Lspe-U0-Lsp0)+pres_init0*A;
normNRJ = Em(1);
diff_M  = max(abs((M-M(1))./M(1)));
fprintf(1,'Fluid mass variation = %g percent \n',diff_M);

titre=['U_o=',num2str(U0),'m - F(0)=',num2str(F0,3),' N - E_m(0)=',num2str(normNRJ,3),' J - K_s=',num2str(vprel(1),3),' N/m,  M_s=',num2str(vprel(2)),' kg'];
indice= ceil(linspace(1,npas,10)); 
%
% displaying pressures
%
figure(4);
hold on;
ta=[0 t];
presL2t_plot=[pres_init.*[1 1 1]; presL2t];
plot(ta./T0,presL2t_plot(:,1)./1E5,'r--',ta./T0,presL2t_plot(:,2)./1E5,'b--',ta./T0,presL2t_plot(:,3)./1E5,'k--','linewidth',1)
hold on
plot(ta(indice)./T0,presL2t_plot((indice),1)./1E5,'r*',ta(indice)./T0,presL2t_plot((indice),2)./1E5,'bs',ta(indice)./T0,presL2t_plot((indice),3)./1E5,'ko','markersize',10)
ylabel('Piston pressure (10^5 Pa)')
xlabel('time (s)')


figure(3);
clf
maxfig=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PISTON MOTION
subplot(maxfig,1,1);
plot(t./T0,histo_deformation(:,1)./U0,'k-','linewidth',1);
axis tight;
ylabel('U(t)/U0');
title(titre,'fontsize',11);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FORCE
subplot(maxfig,1,2)
dQdM_dt([2:npas])=diff(QdM)./Delta_t_storage([2:npas]);dQdM_dt(1)=dQdM_dt(2);
plot(t./T0,-dQdM_dt./F0,'k',t./T0,Force_ext./F0,'b-','linewidth',1)
hold on
plot(t(indice)./T0,-dQdM_dt(indice)./F0,'ko',t(indice)./T0,Force_ext(indice)./F0,'bs','linewidth',1)
axis tight
ylabel('Norm. Forces','fontsize',11)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ENERGIES
subplot(maxfig,1,3)
plot(t./T0,-Imp_fl(:,indice_fluid_model)./normNRJ,'k-','linewidth',1)
hold on
plot(t(indice)./T0,-Imp_fl(indice,indice_fluid_model)./normNRJ,'ko')
%
plot(t(indice)./T0,(Em(indice)-0*Em(1))./normNRJ,'bs-','linewidth',1)
%
plot(t./T0,(Ec)./normNRJ,'r--')%,t./T0,(Ep-0*Em(1))./normNRJ,'m-.')
axis tight
ylabel('Norm. Energies','fontsize',11)
xlabel('t/T_o','fontsize',10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Max_Em=abs(Em(npas)-Em(1));
Max_Im=abs(Imp_fl(npas,indice_fluid_model));
diff_NRJ=abs(Em(npas)-Em(1)-Imp_fl(npas,indice_fluid_model))./abs(Imp_fl(npas,indice_fluid_model));
fprintf(1,'Energy gap between fluid and structure = %g percent \n',diff_NRJ);

if ((fluid_model == 'C')&(i3Dview==1))
    figure(1);
    plot3(courbet,courbex,pres_init(1)./1e5.*ones(length(courbex),1),'b-');
    axis tight;
    xlabel('Time (s)','fontsize',12);
    ylabel('x (m)','fontsize',12);
    zlabel('Chamber Pressure (10^5 Pa)','fontsize',12);
    view([-132, 40]);
    grid;
end
