% CALCULATE RELATIVE POSITION DATA (from excitation to screen)

function [sc,r] = generateRelPositionData(N,b,lambda,twopi);

% Establish sine and cosine values for vectors from one 
% position in x0 to all positions in x1, and find distances 
% and relative phase differences between the points. 
% Function is written by AIV. Version 15. October 2017

y = [-N:N]; 
b2 = b*b*1.0;  
y2p = (y.*y) + b2;
rnn = sqrt(y2p);
sc(:,1) = b./rnn;
sc(:,2) = y./rnn;
r(:,1) = 1./sqrt(rnn);
fs = mod(rnn,lambda);
r(:,2) = fs.*(twopi/lambda);
% mx = max(r(:,1));  % For testing if field reduction vs 
					 % distance is correct
% r(:,1) = mx;
%  plot(x2,r(:,2),'-k'); % Test plot of these variables 
%  figure;
return;