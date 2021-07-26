% MAKE AN ANIMATION OF ALL SPATIAL WAVES}

function [xavt] = ...
   pg_animer(x,deltat,N,A,phase,k,omega,imin,imax,xmax,disp)

% Animation of a wave package during some time. To ease the 
% understanding of the difference between phase and group 
% velocity, a pure monochromatic wave with the central 
% wavelength is animated along with the wave package.
% Returns how far the monochromatic wave has moved during 
% the animation (indicates the phase velocity).
% Input parameters: See the explanations given in the
% functions pg3.m, pd:wpack.m, pg_fft.m, pg_omega.m and 
% pg:wave.m  Version: Oct 5 2017 AIV

figure;
count=1;
% The animation loop
for n = 1:200
    % Calculate the wave at time t (manual IFFT)
    t = deltat*n;
    [zrecon] = pg_wave(x,t,N,A,phase,k,omega,imin,imax);
    % Calculate also the wave with central spatial frequency 
    % in the distribution
    imean = round((imin+imax)/2.0);
    [zrecon0] = pg_wave(x,t,N,A,phase,k,omega,imean,imean);
    % Calculate marking positions, start and end of movement 
    % at phase velocity
    x00 = xmax/8.0;
    xavt = x00 + t*omega(imean)/k(imean);
    % Plots everything
    plot(x,2.5*zrecon0,'-g', x,zrecon,'-b', x00,0.25,'+r', ...
       xavt,0.25,'+r');
          xlabel('Position (rel)');
          ylabel('Amplitude (rel)');
          axis([0,xmax,-1.04,1.04])
          title('Movement to a blue wave package');
          S = sprintf('Time: %.2f s',t);
          text(3.0, 0.8,S);
          S = sprintf('Xref: %.2f',xavt);
          text(3.0, 0.65,S);
          S = sprintf('Dispersion code: %.1f',disp);
          text(3.0, -0.8,S);
          M(count)=getframe;
          count=count+1;
          M(count)=getframe;
          count=count+1;
end
% Animation is played with (1 x 20 frames per sec)
movie(M,1,20);