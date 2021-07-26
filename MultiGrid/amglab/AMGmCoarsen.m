function [toodeep] = AMGmCoarsen(level)
% toodeep will be set to 1 if the next coarse level could not be set up, otherwise

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

% Uses subroutines written by Gundolf Haase

amg_globals;

if DEBUG == 1
    disp('Starting AMGmCoarsen')
end

% if the edge and edge connectivity matrices are empty for this level
if isempty( A(level).edges ) % determine them...
    % derive node numbering for edge matrices from element connectivity
    [Ee,ea] = GetEdgeElements_Schur( A(level).elements , A(level).elconn , PAIRS);
    % Ee is edge matrix having dimensions [2,2,n_edges]
    % ea is edge connectivity matrix having dimensions [2,n_edges]
else % bring them in for use in this algorithm
    Ee = A(level).edges;
    ea = A(level).edconn;
end

%  Determine neighborhood for nodes which are connected via an edge.
%  neigh(i).nn is the array of neigbors for node i (incl. node i).
neigh = NodeNeighborsEdges(ea);

% ----------------------- INITIAL COARSENING PHASE -----------------------% 

% initialize coarse, working, and fine sets to be empty
C_set = [];
F_set = [];

nnodes = size(neigh); % get the number of nodes in the fine grid

% initialize the working set to all fine grid nodes
U_set=[];
for m=1:nnodes(2) % for each node in the fine grid
    U_set = union(U_set, m); % add this node to the working set
end

for m=1:nnodes(2) % for each node in the fine grid
    SC(m).set = strong(m, U_set, level); % get the set of strongly connected nodes to node m
    SC(m).phi = size( SC(m).set ); % get the cardinality of the strongly connected set of node m
end

n=0; % counts the number of nodes which have not been marked coarse or fine
loopcount = 0;
while n < nnodes(2) % while there are nodes that have not been marked
    loopcount = loopcount+1;
    maxPHI = 0; i=0; % initialize node with maximal strongly connected set
    usize = size( U_set ); % get the size of the working set
    
    % find the node i in the working set with maximal strongly connected set
    for m=1:usize(2) % for each node m in the working set
        if SC( U_set(m) ).phi > maxPHI % if the strogly connected set of m is the maximum encountered
            i=U_set(m); % mark the node having the max value
            maxPHI = SC( U_set(m) ).phi; % remember the max value
        end
    end
    
    C_set = union(C_set, i); % add the node with maximal strongly connected set to the coarse grid
    U_set = setdiff(U_set, i); % remove that node from the working set
    n = n+1; % count that one more node has been marked
    
    if i ~= 0 % if a suitable node i is found
       if SC(i).phi ~= 0 % if the set of strongly connected points to i is non-empty
          SiU = intersect( SC(i).set , U_set ); % construct the set of strongly connected nodes to i which are still in the working set
       end
    else % if a suitable node i is not found
          SiU = []; % the intersect is empty
    end
    
    ssize = size( SiU ); % get the size of that set
    for j=1:ssize(2) % for each node in that set
        F_set = union(F_set, SiU( j ) ); % tentatively mark that node as fine
        U_set = setdiff(U_set, SiU( j ) ); % remove that node from the working set
        n = n+1; % count that one more node has been marked
        
        SjU = intersect( SC(j).set , U_set ); % construct the set of strongly connected nodes to j which are still in the working set
        sssize = size( SjU ); % get the size of that set
        for k=1:sssize(2) % for each node in that set
            SC( SjU(k) ).phi = SC( SjU(k) ).phi + 1; % incriment the number of strong connection to that node
        end
    end 
end
    
% ------------------------ FINAL COARSENING PHASE ------------------------% 

% skipped in this version...

% ------------------- GENERATE INTERPOLATION MOLECULES -------------------% 

fnodes = size(F_set); % get the number of selected fine nodes
cnodes = size(C_set); % get the number of selected coarse nodes

% set up storage for global interpolation operator
% zero matrix with dimensions (all nodes on fine grid) X (coarse selected nodes)
W(level).Iweight = zeros( nnodes(2), cnodes(2) );

% generate interpolation molecules for each node selected as fine (row in molecule matrix)
for ni=1:fnodes(2) % for each marked fine node ni
    NNsize = size( neigh( F_set(ni) ).nn, 2 ); % get the size of the node neighborhood around the i-th fine point
    NEsize = size( neigh( F_set(ni) ).ne, 2 ); % get the size of the edge neighborhood around the i-th fine point

    % make sets of C and F nodes in the neighborhood
    fineneighbors = intersect( neigh( F_set(ni) ).nn, F_set );
    fncount=size(fineneighbors);
    coarseneighbors = intersect( neigh( F_set(ni) ).nn, C_set );
    cncount=size(coarseneighbors);
    
    % generate a local node numbering for the neighborhood of i
    % fine points should be sorted to the front of the array
    loc_glob(1:fncount(2))=fineneighbors;
    loc_glob(fncount(2)+1:NNsize)=coarseneighbors;
    % find the local numbering of node ni
    [is,i] = ismember( F_set(ni), loc_glob );
    
    % initialize molecule components for this fine node to the appropriate sizes
    M( ni ).ff = zeros(fncount(2), fncount(2));
    M( ni ).fc = zeros(fncount(2), cncount(2));
    % M( ni ).cf = zeros(cncount(2), fncount(2));
    % M( ni ).cc = zeros(cncount(2), cncount(2));

    for ecount=1:NEsize % for each edge neighboring node i
        % get the global index in the connectivity matrix for this edge
        e = neigh( F_set(ni) ).ne(ecount);
        
        i1 = ea(1,e); % find the first endpoint of this edge in global numbering
        i2 = ea(2,e); % find the second endpoint of this edge in global numbering
        
        i1C=0;i1F=0;i2C=0;i2F=0; % initialize flags
        % attempt to find i1's local numbering from the coarse node neighbors of i
        [i1C, i1_loc] = ismember(i1, coarseneighbors);
        if ~i1C % if i1 is not a local course point
            % find its local numbering from the fine node neighbors of i
            [i1F, i1_loc] = ismember(i1, fineneighbors);
        end
        
        % attempt to find i2's local numbering from the coarse node neighbors of i
        [i2C, i2_loc] = ismember(i2, coarseneighbors);
        if ~i2C % if i1 is not a local course point
            % find its local numbering from the fine node neighbors of i
            [i2F, i2_loc] = ismember(i2, fineneighbors);
        end
        
        % place the edge matrix components of this edge into the appropriate molecule component
        if [i1F i2F] % if both endpoints are finepoints
            % add the edge matrix contributions into the Mcc molecule component
            M( ni ).ff(i1_loc,i1_loc) = M( ni ).ff(i1_loc,i1_loc) + Ee(1,1,e);
            M( ni ).ff(i1_loc,i2_loc) = M( ni ).ff(i1_loc,i2_loc) + Ee(1,2,e);
            M( ni ).ff(i2_loc,i1_loc) = M( ni ).ff(i2_loc,i1_loc) + Ee(2,1,e);
            M( ni ).ff(i2_loc,i2_loc) = M( ni ).ff(i2_loc,i2_loc) + Ee(2,2,e);
        % elseif [i1C i2C] % if both endpoints are coarse points
        %     % add the edge matrix contributions into the Mcc molecule component
        %     M( ni ).cc(i1_loc, i1_loc) = M( ni ).cc(i1_loc, i1_loc) + Ee(1,1,e);
        %     M( ni ).cc(i1_loc, i2_loc) = M( ni ).cc(i1_loc, i2_loc) + Ee(1,2,e);
        %     M( ni ).cc(i2_loc, i1_loc) = M( ni ).cc(i2_loc, i1_loc) + Ee(2,1,e);
        %     M( ni ).cc(i2_loc, i2_loc) = M( ni ).cc(i2_loc, i2_loc) + Ee(2,2,e);
        elseif i1C % endpoint one is coarse and two is fine
            % add the edge matrix contributions into the molecule
            % M( ni ).cc(i1_loc, i1_loc) = M( ni ).cc(i1_loc, i1_loc) + Ee(1,1,e);
            % M( ni ).cf(i1_loc, i2_loc) = M( ni ).cf(i1_loc, i2_loc) + Ee(1,2,e);
            M( ni ).fc(i2_loc, i1_loc) = M( ni ).fc(i2_loc, i1_loc) + Ee(2,1,e);
            M( ni ).ff(i2_loc, i2_loc) = M( ni ).ff(i2_loc, i2_loc) + Ee(2,2,e);
        else % endpoint one is fine and two is coarse
            % add the edge matrix contributions into the molecule
            M( ni ).ff(i1_loc, i1_loc) = M( ni ).ff(i1_loc, i1_loc) + Ee(1,1,e);
            M( ni ).fc(i1_loc, i2_loc) = M( ni ).fc(i1_loc, i2_loc) + Ee(1,2,e);
            % M( ni ).cf(i2_loc, i1_loc) = M( ni ).cf(i2_loc, i1_loc) + Ee(2,1,e);
            % M( ni ).cc(i2_loc, i2_loc) = M( ni ).cc(i2_loc, i2_loc) + Ee(2,2,e);
        end
        
        % see if there are any contributions from coarse edges not connected directly to i
        
        % generate the set of fine node neigbors to i, excluding i
        FNI = setdiff( fineneighbors, i );
        FNIsize = size( FNI, 2 );% get the size of that set
        for j=1:FNIsize % for each non-i fine direct neighbor of i, called j
            [is,j_loc] = ismember( FNI( j ), fineneighbors ); % local numbering of j
            for k=1:cncount(2) % and for each coarse neighbor of i, called k
                k_loc=k; % local numbering of k
                % find the edges which connect the fine direct neighbor 
                % of i, called j, to the coarse direct neighbor of i, called k
                % NOTE: such an edge is not guaranteed to exist, there can be at most one
                IND = intersect( neigh( FNI(j) ).ne , neigh( coarseneighbors(k) ).ne );
                
                if ~isempty(IND) % if such an edge was found between j and k
                    ne = IND(1); % get the global indedx for the connectivity matrix of the edge
                    
                    % set up the indices for the edge matrices
                    if ea(1,ne) == loc_glob(k_loc) % if the first endpoint is our coarse point k
                    % NOTE: in this case the second endpoint must be our fine node j
                        ee1=1; ee2=2;
                    else % endpoint one is fine and two is coarse
                        ee1=2; ee2=1;
                    end
                    % add the edge matrix contributions into the molecule
                    % M( ni ).cc(k_loc, k_loc) = M( ni ).cc(k_loc, k_loc) + Ee(ee1,ee1,ne);
                    % M( ni ).cf(k_loc, j_loc) = M( ni ).cf(k_loc, j_loc) + Ee(ee1,ee2,ne);
                    M( ni ).fc(j_loc, k_loc) = M( ni ).fc(j_loc, k_loc) + Ee(ee2,ee1,ne);
                    M( ni ).ff(j_loc, j_loc) = M( ni ).ff(j_loc, j_loc) + Ee(ee2,ee2,ne);
                end
            end
        end 
    end
 
    % store the local to global numbering utilized in this ni-th molecule for later use
    % M( ni ).numbering = loc_glob;

    % generate the fncount by cncount matrix Pfc which will be used as the
    % interpolation weight matrix for this fine point
    M( ni ).Iweight = (-(inv( M( ni ).ff ))) * M( ni ).fc;
    
    % extract the i_th row of Pfc to go into the global interpolation operator
    for cc=1:cncount(2)
        globalcoarse = loc_glob(fncount(2)+cc);
        [is,ci] = ismember( globalcoarse, C_set );
        W(level).Iweight( F_set(ni), ci ) = M(ni).Iweight(i,cc);
    end
end

% -------------------------- COARSE LEVEL SETUP --------------------------% 

% generate coarse edges
currentCedge = 0;

% get the number of coarsely selected points
cnodes = size(C_set);

% for each selected coarse point i
for icount=1:cnodes(2)
    i = C_set( icount ); % get the global point ID for i
    
    % for each selected coarse point j
    for jcount=icount+1:cnodes(2)
        j = C_set( jcount ); % get the global point ID for j
        
        % examine the possibility of making a coarse edge between i and j
        if i ~= j % do not consider cases where it is a self connection...
            
            CE(icount,jcount).cc = zeros(2,2);
            CE(icount,jcount).cf = zeros(2,2);
            CE(icount,jcount).fc = zeros(2,2);
            
            % if i and j are directly connected on the fine grid
            if ismember( j, neigh(i).nn )
                % create a coarse edge between i and j
                currentCedge = currentCedge + 1;
                % place the edge in the connectivity matrix
                A(level+1).edconn(1,currentCedge) = icount; % coarse grid numbering
                A(level+1).edconn(2,currentCedge) = jcount; % coarse grid numbering
                
                % find the edge connecting i and j on the fine grid
                end_i=-1; end_j=-1;
                for( edgenum=1:size(ea,2) ) % check all fine endges
                    [is loc] = ismember([i;j], ea(:,edgenum), 'rows');
                    if is % if current fine edge is ij or ji
                        if loc(1)==1 % if current fine edge is ij
                            end_i = 1; % ea(1,edgenum);
                            end_j = 2; % ea(2,edgenum);
                        else % current fine edge is ji
                            end_i = 2; % ea(2,edgenum);
                            end_j = 1; % ea(1,edgenum);
                        end
                        break; % stop looking for the edge
                    end
                end
                
                % add the contribution of that edge to the coarse edge molecule
                % sprintf('end_i: % d end_j: % d',end_i,end_j)
                CE(icount,jcount).cc(1,1) = Ee(end_i, end_i, edgenum);
                CE(icount,jcount).cc(1,2) = Ee(end_i, end_j, edgenum);
                CE(icount,jcount).cc(2,1) = Ee(end_j, end_i, edgenum);
                CE(icount,jcount).cc(2,2) = Ee(end_j, end_j, edgenum);
            else % i and j are not directly connected on the fine grid            
                % for each fine selected point k
                fnodes = size(F_set);
                for kcount=1:fnodes(2)
                    k = F_set( kcount ); % get the global point ID for k
                
                    % construct the set of coarse points which are strongly connected to fine point k
                    Sck = intersect( C_set, SC(k).set );
                
                    % if k is connected to both coarse points i and j
                    if ismember( [i j] , Sck )
                        % create a coarse edge between i and j
                        currentCedge = currentCedge + 1;
                        % place the edge in the connectivity matrix
                        A(level+1).edconn(1,currentCedge) = icount; % coarse grid numbering
                        A(level+1).edconn(2,currentCedge) = jcount; % coarse grid numbering
                        
                        % find the edge connecting i and k on the fine grid
                        end_i=-1; end_k=-1;
                        for( edgenum=1:size(ea,2) ) % check all fine endges
                            [is loc] = ismember([i;k], ea(:,edgenum), 'rows');
                            if is % if current fine edge is ik or ki
                                if loc(1)==1 % if current fine edge is ik
                                    end_i = 1; % ea(1,edgenum);
                                    end_k = 2; % ea(2,edgenum);
                                else % current fine edge is ki
                                    end_i = 2; % ea(2,edgenum);
                                    end_k = 1; % ea(1,edgenum);
                                end
                                break; % stop looking for the edge
                            end
                        end                        
                        
                        % add the contribution of that edge to the coarse edge molecule
                        if end_i > 0 && end_k > 0
                        CE(icount,jcount).cc(1,1) = Ee(end_i, end_i, edgenum);
                        CE(icount,jcount).cf(1,2) = Ee(end_i, end_k, edgenum);
                        CE(icount,jcount).fc(2,1) = Ee(end_k, end_i, edgenum);
                        CE(icount,jcount).ff(2,2) = Ee(end_k, end_k, edgenum);
                        end
                        % find the edge connecting k and j on the fine grid
                        end_k=-1; end_j=-1;
                        for( edgenum=1:size(ea,2) ) % check all fine endges
                            [is loc] = ismember([k;j], ea(:,edgenum), 'rows');
                            if is % if current fine edge is kj or jk
                                if loc(1)==1 % if current fine edge is kj
                                    end_k = 1; % ea(1,edgenum);
                                    end_j = 2; % ea(2,edgenum);
                                else % current fine edge is jk
                                    end_k = 2; % ea(2,edgenum);
                                    end_j = 1; % ea(1,edgenum);
                                end
                                break; % stop looking for the edge
                            end
                        end                     
                        
                        % add the contribution of that edge to the coarse edge molecule
                        if end_k > 0 && end_j > 0
                        CE(icount,jcount).ff(1,1) = Ee(end_k, end_k, edgenum);
                        CE(icount,jcount).fc(1,2) = Ee(end_k, end_j, edgenum);
                        CE(icount,jcount).cf(2,1) = Ee(end_j, end_k, edgenum);
                        CE(icount,jcount).cc(2,2) = Ee(end_j, end_j, edgenum);
                        end
                    end
                end
            end                 
        end
    end
end


% determine the number of course edges which were selected
numcourseedges = size(A(level+1).edconn,2);
if( numcourseedges < 3 ) % not enough course edges are selected to form an element
% NOTE: it would also be good to make sure that the 3 or more edges actually
% form an element. I propose using a simple Kruskhal's algorithm to
% determine if a circular link is present between the existing edges
    toodeep = 1;
else
    toodeep = 0;
end

% generate the coarse edge matrices
for( CEcount=1:size(A(level+1).edconn,2) )
    i=A(level+1).edconn(1,CEcount);
    j=A(level+1).edconn(2,CEcount);
    A(level+1).edges(:,:,CEcount) = CE(i,j).cc - ( CE(i,j).cf * CE(i,j).fc );
end

% generate coarse grid matrix via the Edge Accumulation Method
A(level+1).matrix = AccumulateMatrix( A(level+1).edges , A(level+1).edconn );
tester2 = A(level+1).matrix;

% generate coarse grid matrix via the Galerkin approach
W(level).Rweight = W(level).Iweight'; % restriction is the transpose of the interpolation operator in this case
A(level+1).matrix = W(level).Rweight * A(level).matrix * W(level).Iweight;
tester1 = A(level+1).matrix;

if DEBUG == 1
    disp('Exit AMGmCoarsen')
end

% difference = abs( tester1 - tester2 )