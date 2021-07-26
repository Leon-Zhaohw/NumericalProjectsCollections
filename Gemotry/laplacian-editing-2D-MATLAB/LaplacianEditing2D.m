%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script for 2D Laplacian Editing
%
% Please refer to "Laplacian Surface Editing" by Olga Sorkine,
%   Daniel Cohen-Or, Yaron Lipman,  Marc Alexa, Christian Rössl and
%   Hans-Peter Seidel,
%   Eurographics/ACM SIGGRAPH Symposium on Geometry Processing 2004,
%   pp. 179--188, ACM Press.
%   http://www.cs.nyu.edu/~sorkine/ProjectPages/Editing/lse.html
%
% The script demonstrates how to build the Laplacian Editing matrix
% when editing 2D curves. The static and handle vertices are hard-coded
% for some example curves, provided with this script.
% The script prompts the user to input new locations for the handle
% vertices and displays the result in a new figure.
%
% Note: only the first stage of Laplacian editing is computed, namely
% editing with implicit local similarity transformations. The second stage
% (renormalizing the delta-coordinates to remove the introduced local
% scaling) is not implemented here.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assume xy is n x 2

% uncomment the curve you would like to try
%xy = load('unicorn.txt');
xy = load('T.txt');
%xy = load('sine_circle.txt');


n=length(xy);

% the Laplacian matrix (uniform weighting)
L = spdiags(ones(n,1),0,n,n) - spdiags(ones(n,1),1,n,n);
L = L+L';
L(1,n)= -1;
L(n,1) = -1;
L = L./2;


delta = L*xy;


% we want to construct the matrix of the system for v-primes
L_prime = [   L     zeros(n)   % the x-part
	       zeros(n)    L    ]; % the y-part
	      
for i=1:n
  % the neighbors of i are i-1 and i+1
  ring = [(1 + mod(i-2,n)), i, (1 + mod(i,n))]; % i-1, i, i+1
  V = xy(ring,:)';
  V = [V
       ones(1,length(ring))];


  % the coeff matrix for the system that solves for T
  %      s  a  t1
  % T = -a  s  t2
  %      0  0  1
  
  C = zeros(6,4);
  % ... Fill C in
  for r=1:length(ring)
    C(r,:) =                [V(1,r)       V(2,r)  V(3,r)      0  ];
    C(length(ring)+r,:) =   [V(2,r)  (-1)*V(1,r)       0  V(3,r) ];

  end;
    
  Cinv = pinv(C);
  s =  Cinv(1,:);
  a =  Cinv(2,:);

  delta_i = delta(i,:)';
  delta_ix = delta_i(1);
  delta_iy = delta_i(2);
  
  % T*delta gives us an array of coefficients
  Tdelta = [delta_ix*s      + delta_iy*a 
	        delta_ix*(-1)*a + delta_iy*s];
	    
  
  % updating the weights in Lx_prime, Ly_prime, Lw_prime
  L_prime(i,[ring (ring + n)]) = L_prime(i,[ring (ring + n)]) +...
                                              (-1)*Tdelta(1,:);
  L_prime(i+n,[ring (ring + n)]) = L_prime(i+n,[ring (ring + n)]) +...
                                                (-1)*Tdelta(2,:);
end;


% additional constraints - we need at least 3

% weight for the constraints
w=1;

% in the following, uncomment the appropriate curve 

%static_anchors = [196:198 114:116]; % unicorn
%handle_anchors = [152 136];         % unicorn

static_anchors = [157:159 470:472]; % sine_circle
handle_anchors = [314];             % sine_circle

%static_anchors = [99:101 240:242]; % T
%handle_anchors = [171];            % T

%static_anchors = [99:101 240:242];  % T with two handle anchors
%handle_anchors = [171 172];         % T with two handle anchors


% building the least-squares system matrix
A_prime = L_prime;
rhs = zeros(2*n,1);
anch_pos = [];

anchors = [static_anchors handle_anchors];
for j=1:length(anchors)
  A_prime = [A_prime
	     w*((1:(2*n))==anchors(j))
	     w*((1:(2*n))==(anchors(j)+n))];
  rhs = [rhs
	 w*xy(anchors(j),1)
	 w*xy(anchors(j),2)];
  
  anch_pos = [anch_pos
	      xy(anchors(j),1:2)];
end;


% displaying the curve: static anchors are in black, the handle to be moved
% in red.
figure(1);
plot(xy(:,1),xy(:,2),'b-',...
     xy(static_anchors,1),xy(static_anchors,2),'*k',...
     xy(handle_anchors,1),xy(handle_anchors,2),'*r');
axis equal;

% moving the handle (ok to click outside the figure axes as well)
[x_input y_input] = ginput(1);
lenr = length(rhs);
rhs((lenr-1):lenr) = [x_input y_input]';
anch_pos(length(anch_pos),:) = [x_input y_input];

% moving second handle vertex if exists
if length(handle_anchors) > 1
    [x_input y_input] = ginput(1);
    rhs((lenr-3):(lenr-2)) = [x_input y_input]';
    anch_pos(length(anch_pos)-1,:) = [x_input y_input];
end


% solving for v-primes
xy_col = A_prime\rhs;
xy_prime_nonh = [xy_col(1:n) xy_col((n+1):(2*n))];

figure(2);
plot(xy(:,1),xy(:,2),'g-',...
     xy_prime_nonh(:,1),xy_prime_nonh(:,2),'-b',...
     xy_prime_nonh(anchors,1),xy_prime_nonh(anchors,2),'*k',...
     anch_pos(:,1),anch_pos(:,2), 'mo',...
     xy_prime_nonh(anchors(length(anchors)),1),xy_prime_nonh(anchors(length(anchors)),2),'*r');
axis equal;
title('Editing result');


