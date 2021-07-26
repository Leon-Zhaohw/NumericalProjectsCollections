% This is a scrip to inforce the Dirichlet boundry condition on FEM
% problems.  At the moment this script will only work on a 50 element
% system.  Once I am sure that it works, I will make it more general.


amg_globals;

[a b] = size(A(1).matrix);

a = (sqrt(a)-2)^2;

M1 = zeros( a , b );

k = sqrt(b) + 2;

for i = 1:a
    M1(i,:) = A(1).matrix(k, :);
    k = k + 1;
    if mod(k,sqrt(b)) == 0
        k = k + 2;
    end
end

M2 = zeros(a);

k = sqrt(b) + 2;

for i = 1:a
    M2(:,i) = M1(:,k);
    k = k + 1;
    if mod(k,sqrt(b)) == 0
        k = k + 2;
    end
end

A(1).matrix = M2;

clear M1;
clear M2;

A(1).elements = [];
A(1).elconn = [];

dim = sqrt(a);
    nodecount = 0;
    % for each sqare in the stencil (two linear triangle elements)
    for i=1:dim-1
        for j=1:dim-1
            idx = (i-1)*dim + j;
            nodecount=nodecount+1; % generate a new element
            % generate and store edge matrix for the first element in this square
            A(1).elements(:,:,nodecount) = [1,0,-1;0, 1, -1; -1, -1, 2];
            % generate and store connectivity matrix for the first element in this square
            A(1).elconn(:,nodecount) = [idx,idx+dim+1,idx+dim];
            
            nodecount=nodecount+1; % generate a new element
            % generate and store edge matrix for the second element in this square
            A(1).elements(:,:,nodecount) = [1, -1, 0; -1, 2, -1; 0, -1, 1];
            % generate and store connectivity matrix for the second element in this square
            A(1).elconn(:,nodecount) = [idx,idx+1,idx+dim+1];
            %disp('A(1).elements in loop')
            %A(1).elements
            %disp('A(1).elconn in loop')
            %A(1).elconn
        end
    end

    %A(1).elements
    
%A(1).elconn = [1 1 2 2 3 3  5  5  6  6  7  7  9  9 10 10 11 11;
%               6 2 7 3 8 4 10  6 11  7 12  8 14 10 15 11 16 12;
%               5 6 6 7 7 8  9 10 10 11 11 12 13 14 14 15 15 16];
           