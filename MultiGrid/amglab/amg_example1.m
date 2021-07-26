% AMG_EXAMPLE1 sets us up for Poisson Equation

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

function amg_example1;

amg_globals;

if(PROB_TYPE==STIFFNESS) % if we are setting up is a finite difference stiffness matrix
    % FINEPOINTS is the dimension of the square matrix
    dim = sqrt(FINEPOINTS); % dim is the square root of that so we can use gallery

    % construct the stiffness matrix
    % a dim^2 by dim^2 (FINEPOINTS x FINEPOINTS) possion discretization
    % 9-point stencil used
    A(1).matrix = full( gallery('poisson', dim));
elseif(PROB_TYPE==ELEMENTS)
    MAXVERTEX=3; % linear triangle
    PAIRS = int16([ 1 2; 1 3; 2 3 ]);%  edges in lin. triangle
    
    % initialize edge matrices
    A(1).elements = zeros(MAXVERTEX, MAXVERTEX, FINEPOINTS);
    %disp('A(1).elements before loop')
    %A(1).elements
    % initialize element connectivity matrix
    A(1).elconn = zeros(MAXVERTEX, FINEPOINTS);
    %disp('A(1).elconn before loop')
    %A(1).elconn
    
    dim = sqrt(FINEPOINTS/2) + 1;
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
    
    % initialize the edge and edge connectivity matrices
    A(1).edges = [];
    A(1).edconn = [];

    % accumulate stiffness matrix from the element matrices
    A(1).matrix = AccumulateMatrix( A(1).elements , A(1).elconn );
    %A(1).matrix
    %size(A(1).matrix)
    
    % tester=A(1).matrix
    % tester2=full( gallery('poisson', 3))

    Dirichlet_Correction;
    
end



dim = size(A(1).matrix);
    
X_Guess = zeros( dim(1), 1 );

RHS = ones( dim(1), 1 );