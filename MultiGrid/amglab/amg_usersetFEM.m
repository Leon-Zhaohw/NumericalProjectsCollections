%AMG_USERSET_FEM sets us up for User Specified Example using Finite Element
%Discretization

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

%Adapted from a code by Gundolf Haase
%Algorithms provided by Johannes Kraus

amg_globals;

filename = input('\nPlease enter the file name containing your FEM matrix.\n(Press Return for the Default):', 's');
if isempty(filename)
    filename = 'e32_isb.txt';
end

fin = fopen (filename,'rt');
FINEPOINTS = fscanf(fin,'%4d',1);         % number of elements
MAXVERTEX = fscanf(fin,'%4d',1);      % max. number of vertices per element
if MAXVERTEX==3 % linear triangle
    PAIRS = int16([ 1 2; 1 3; 2 3 ]);% edges in lin. triangle
else % linear tetrahedron (3 and 6 are only possibilities)
    PAIRS = int16([ 1 2; 1 3; 2 3; 1 4; 2 4; 3 4 ]);
end;
%maxedge = length(PAIRS);

%initialize edge matrices
A(1).elements = zeros(MAXVERTEX, MAXVERTEX, FINEPOINTS);
%initialize element connectivity matrix
A(1).elconn = zeros(MAXVERTEX, FINEPOINTS);
for i=1:FINEPOINTS
    fscanf(fin,'%s',1); %read in node header and discard
    ii = fscanf(fin,'%d',1); fscanf(fin,'%1s',1); %read in node number
    if i~=ii,  exit; end; %compare node number to current node (file integrity)
    %read in edge matrix for this node
    A(1).elements(:,:,i) = fscanf(fin,'%10g',[MAXVERTEX, MAXVERTEX]);
    %read in connectivity matrix for this node
    A(1).elconn(:,i)   = fscanf(fin,'%4d',MAXVERTEX);
end;
fclose(fin);

%initialize the edge and edge connectivity matrices
A(1).edges = [];
A(1).edconn = [];

%accumulate stiffness matrix from the element matrices
A(1).matrix = AccumulateMatrix( A(1).elements , A(1).elconn );

dim = size(A(1).matrix);

X_Guess = zeros( dim(2), 1 );

RHS = ones( dim(2), 1 );