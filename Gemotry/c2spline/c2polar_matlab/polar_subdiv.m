% Helper function for c2_polar_spline.m

% Input: 2-layer polar configuration polcfg with 4n vertices,
%        where polcfg(1:n) is the same vertex
% Output: 4-layer polar configuration q with 6n vertices.
% From Myles, Karciauskas, Peters.
%      Extending Catmull-Clark Subdivision and PCCM with Polar Structures.
function q = polar_subdiv(n, polcfg)

beta = 5.0/8.0;
alpha = beta - 1.0/4.0;
c = cos(2*pi/n * [0:n-1]);
gamma = 1/n * (beta - 1.0/2.0 + 5.0/8.0*c + c.^2 + 1.0/2.0*c.^3);
GammaMatrix = toeplitz(gamma([1,n:-1:2]), gamma);

q = zeros(6,n);
q(1,:) = (1-alpha)*polcfg(1,:) + alpha/n*sum(polcfg(2,:));         % 0-link
q(2,:) = (1-beta)*polcfg(1,:) + polcfg(2,:)*GammaMatrix';          % 1-link
q(3,:) = 0.125*polcfg(1,:) + 0.75*polcfg(2,:) + 0.125*polcfg(3,:); % 2-link
q(4,:) = 0.5  *polcfg(2,:) + 0.5 *polcfg(3,:);                     % 3-link
q(5,:) = 0.125*polcfg(2,:) + 0.75*polcfg(3,:) + 0.125*polcfg(4,:); % 4-link
q(6,:) = 0.5  *polcfg(3,:) + 0.5 *polcfg(4,:);                     % 5-link

