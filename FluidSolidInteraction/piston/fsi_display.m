% ----------------------
% filename: fsi_display.m
%
% The script managing plotting
%
figure(2);
hold off;
subplot(3,1,1);
x=vcor_np1(:,1);
plot(x,vpres./1E5,'b-o');
hold on;
plot(x(nnt),presL2t(istep,1)./1E5,'gs');
plot(x(nnt),presL2t(istep,2)./1E5,'rd');
ylabel('Pressure(++)');
hold off;
titre=['Time : ',num2str(Total_time),' sec.'];
title(titre);

subplot(3,1,2);
plot(x,wx,'--',x,vsol(:,2)./vsol(:,1));
ylabel('w_x, u');

subplot(3,1,3);
plot(x,vsol(:,1));
xlabel('x (m)');
ylabel('\rho');
%
% for displaying shifted profiles as time (3D plottings)
% only valuable for fluid_model == 'C'
%
if ((fluid_model == 'C')&(i3Dview==1))
    figure(1);
    % ---------------
    % We display the initial condition
    if(istep==1)
        plot3(0.*ones(nnt,1),vcor0(:,1),pres_init./1e5.*ones(nnt,1),'r-','linewidth',2);
        hold on;
        plot3(0,vcor0(nnt,1),pres_init./1e5,'ro','markersize',8);
    end
    %
    if(mod(istep-0,modulo)==0)
        its=its+1;
        plot3(Total_time./T0.*ones(nnt,1),vcor_n(:,1),vpres./1e5,'b-','linewidth',1)
        hold on
        plot3(Total_time./T0,vcor_n(nnt,1),presL2t(istep,1)./1e5,'g*','markersize',10)
        plot3(Total_time./T0,vcor_n(nnt,1),presL2t(istep,2)./1e5,'rd','markersize',10)
        plot3(Total_time./T0.*[1 1],vcor_n(nnt,1).*[1 1],[pres_init(1)./1e5 presL2t(istep,3)./1e5],'k--')
        plot3(Total_time./T0,vcor_n(nnt,1),presL2t(istep,3)./1e5,'ko','markersize',8)
        courbet(its)=Total_time./T0;
        courbex(its)=vcor_n(nnt,1);
    end
end
pause(1e-10);
