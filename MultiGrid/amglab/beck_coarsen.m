function [toodeep] = beck_coarsen(level)
% Derrick Cerwinsky
% University of Wyoming
%
% This code was written using the algorithm from the paper
% Graph-Based Algebraic Multugrid for Lagrange-Type Finite
% Elements on Simplicial Meshes by R. Beck.
%
% The only input is the level.  The function will always return a 
% value of zero.

amg_globals;


if issparse(A(level).matrix) == false    % Test if A is sparse, make sparse if not.
    A(level).matrix = sparse(A(level).matrix);
end

[m n] = size(A(level).matrix);  % Number of nodes.  The matrix is square, so we do not keep the m.

m = nnz(A(level).matrix);  % Number of non-zero elements in A.

I = zeros( 1, n );  % Initialize the set of nodes.  0 = undeclared, 1 = coarse, 2 = fine.

c = 0;    % Counts coarse nodes.
f = 0;    % Counts fine nodes.

[row, col, a_ij] = find(A(level).matrix);

cnt = zeros(1,n);  % Count of non-zero elements in each row.
dsp = zeros(1,n);  % Displacement of elements in cnt.


% This counts the number of non-zero elements in each row.

for i = 1:n
    cnt(i) = nnz(A(level).matrix(i,:));
end

% This finds displacement of elements in cnt.
dsp(1) = 1;
for i = 1:n-1
    dsp(i+1) = dsp(i) + cnt(i);
end




% Now we determine which nodes are coarse and which are fine.
for i = 1:n
    if I(i) == 0
        I(i) = 1;
        c = c + 1;
        for j = dsp(i):(dsp(i) + cnt(i) - 1)
            if I(row(j)) == 0
                I(row(j)) = 2;
                f = f + 1;
            end
        end
    end
end


% Now we will make the prolongation matrix

P = sparse(zeros(c + f, c));

k = 1;  % k is a dummy index to track the new index of the coarse nodes.

for i = 1:n
    l = 0;
    if I(i) == 1
        %Q = zeros(1,c);
        %Q(k) = 1;
        %P(i,:) = Q;
        P(i,k) = 1;
        k = k + 1;
    else
        for j = dsp(i):(dsp(i) + cnt(i) - 1)
	    if I(row(j)) == 1
                l = l + 1;
            end
        end
        Q = ones(1,c) .* 1/l;
        P(i,:) = Q;
    end
end

A(level + 1).matrix = sparse(P' * sparse( A(level).matrix * P ) );
W(level).Iweight = P;
W(level).Rweight = P';



toodeep = 0;
