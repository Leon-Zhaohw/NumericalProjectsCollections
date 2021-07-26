% GENERATES THE VARIOUS "EXCITATIONS" (single or double slit etc)

function [x0] = generateExcitation(code,lambda,a,N,Nhalf, ...
	twopi,x0)

% Generate or read in excitation data. NOTE: There are 
% specific requirements for the various excitations that  
% can only be changed in the code below. 
% Function is written by AIV. Version 15. October 2017

switch code
   case (1)
      disp('Single slit')
      m = a * lambda / 2;   % Slit is a wavelengths wide
      x0(Nhalf-m:Nhalf+m-1,1) = 1.0;
      %x0(:,2)= [1:N].*0.05; % Phases are modifies so that  
                             % it mimics a ray is not coming 
                             % perpendicular towards the slit.
   case 2
      disp('Gaussian excitation')
      % Intensity
      width = 200*lambda/2.0;   
      dummy = ([1:N]-Nhalf)./width;
      dummy = (dummy.*dummy);
      x0(:,1) = exp(-(dummy));
      % Phase
      R = 1000; % Radius of curvature in # wavelengths
      y = [-Nhalf:Nhalf-1];
      R2 = R*R*lambda*lambda*1.0;
      dist = sqrt((y.*y) + R2);
      fs = mod(dist,lambda);
      x0(:,2) = fs.*(twopi/lambda);
      %figure;  % Plot if wanted
      %plot(x,x0(:,2),'-r');
      
   case 3
      disp('Straight edge')
      % Excitation is a straight edge, illuminated part: 3/4
      x0(N/4:N) = 1.0;

   case 4
      disp('Double slit')
      % For the double slit, use sufficient large b in 
      % 'parameters' in order to get the well known result
	  x0 = zeros(N,2);
	  a = 20*4;
	  d = 200*4;
	  kx = d/2 + a/2;
	  ki = d/2 - a/2;
	  x0(Nhalf-kx+1:Nhalf-kx+a,1) = 1.0;
	  x0(Nhalf+ki:Nhalf+ki+a-1,1) = 1.0;

   case 5
      disp('Reads excitation data from file')
      % (often earlier calculated results.)
      filename = input('Give name on file with excitation data: ', 's');
      fid = fopen(filename,'r');
      x0(:,1) = fread(fid,N,'double');  % Need to know # 
      									% elements
      x0(:,2) = -fread(fid,N,'double'); 
      status = fclose(fid);
      % figure;  % Testplot to check if data was read properly
      % plot(x,xx0(:,1),'-g');
      % figure;
      % plot(x,xx0(:,2),'-r');
      % aa= xx0(Nhalf);
      % aa  % Test print for one single chosen point
      
   otherwise
      disp('Use code 1-5, please.')
end
return;