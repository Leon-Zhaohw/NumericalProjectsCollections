% GENERATE THE DISPERSION RELATION omega(k)}

function [omega,deltat] = pg_omega(imin,imax,k,disp)

% Generate the dispersion relation omega(k). 
% Version Oct 5 2017, AIV
% Input parameters: imin, imax: first and last index that 
% will be used in the function that creates the animation, 
% k: the wavenumber array created by the function pg_fft, 
% disp: -1, 0, or +1 represent normal, no and anomalous 
% dispersion. 
% Returns: omega: the dispersion relation omega(k),  
% deltat: a suitable delta_t for the animation in order to 
% get useful animation/plots.

if (disp==-1)   % Normal dispersion (here vg = vp/2)
    deltat = 0.015;
    omegafactor = 44.0;
    for i = imin:imax
        omega(i) = omegafactor*sqrt(k(i));
    end
end

if (disp==0)   % No dispersion (here vf = const)
    deltat = 0.015;
    omegafactor = 9.5;
    for i = imin:imax
        omega(i) = omegafactor*k(i);
    end
end

if (disp==1)   % Anomal dispersion (here vg = 3vp/2)
    deltat = 0.0065;
    omegafactor = 5.5;
    for i = imin:imax
        omega(i) = omegafactor*(k(i)^1.5);
    end
end

figure;
plot(k(imin:imax),omega(imin:imax),'.-b');
xlabel('k (rel. units)');
ylabel('omega (rel. units)');
return;