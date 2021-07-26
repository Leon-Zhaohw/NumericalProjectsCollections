function M = IQMeasure(e,B,q,bv)

% Divya Bansal
% Department of Computational Sciences
% University of Kentucky

p1 = (bv - q)' * e;
p2 = B * e;

% calculate inner product for numerator
c = 0;                      % intialize the variable c
n = length(p1);          	% get the lenght of the vector p1
for k=1:n                   % start the loop
	c = c + p1(k) * p1(k);	% update c by the k-th product in inner product
end                         % end loop

% calculate inner product for denominator
d = 0;
m = length (p2);
for l = 1:m
    d = d + (e(l,:) * p2');	
end


M = c/d;