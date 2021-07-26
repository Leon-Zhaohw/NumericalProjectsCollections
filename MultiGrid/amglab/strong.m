function Si = strong(i,V,level)
% Determine the set of strong connections
% strong(I,V,level) returns the (sorted) set of strong connections
% of node I from the remaining set of free nodes V 
% with respect to matrix A at this level and the global treshold SW_BOUND.

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

% Adapted from a code by Dr. Gundolf Haase

amg_globals;

if DEBUG == 1
    disp('Starting strong')
end

Ai = abs( A(level).matrix(i,:) );   %  absolute values of row I in coefficient matrix
dim = size(Ai, 2); %  number of coefficients to consider
aa=0; if i > 1,      aa=max(Ai(1:i-1)); end; % get point with maximum coefficient left of diagonal
bb=0; if i < dim,    bb=max(Ai(i+1:dim)); end; % get point with maximum coefficient right of diagonal
aimax = max(aa,bb); % get point with maximum coefficient

sk_set=[]; % initialize the strongly connected set to be empty
for k=V % consider each point in the domain
    if k~=i % do not consider connections to one's self...
       if aimax~=0 % do not consider zero-connectivity
          if sum(Ai(k))/aimax >= SW_BOUND % If the connection between points I and K is "strong enough"
             sk_set=union(sk_set,k); % add K to the strongly connected set of I
          end
       end
    end
end
Si = sort(sk_set); % sort and return the strongly connected set of I

if DEBUG == 1
    disp('Exit strong')
end