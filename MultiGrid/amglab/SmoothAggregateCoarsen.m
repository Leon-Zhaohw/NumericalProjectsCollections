function [toodeep] = SmoothAggregateCoarsen(level)
% toodeep will be set to 1 if the next coarse level could not be set up, otherwise

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

amg_globals;

% ------------------------- INITIALIZATION PHASE -------------------------% 

nnodes = size(A(level).matrix , 1); % get the size of the stiffness matrix

% generate the R_set, all points in the fine domain which are not isolated
R_set = []; % initialize R_set to be empty
for i=1:nnodes % for each point represented in the stiffness matrix
    for j=i:nnodes % for each connection represented in the stiffness matrix
        % NOTE: assumes symmetric stiffness matrix!
        if A(level).matrix(i,j) ~= 0 % if there existst a connection between i and j
            R_set = union( R_set , j ); % add point j into the R_set
        end
    end
end

Original_Rset = R_set; % make a permanent (within function) copy of the R_set initially
neigh_count = 0; % set the number of aggregated neighborhoods to zero
C_set = []; % initialize the C_set to be empty

% ------------------- INITIAL COARSENING (AGGREGATION) -------------------% 

% here we will select disjoint strongly coupled neighborhoods
Static_Rset = R_set; % make a static copy of the R_set initially
Rsize = size(R_set); % find the initial size of the R_set
% divya
% level

i=1; % initialize the point consideration counter

while i <= Rsize(2) % while there are still points to be considered for neighborhood formation (in the static set)
    [dummy,location]=ismember( Static_Rset(i), R_set );
    if dummy % if the point we consider in the static set is still a member of the R_set
        S(Static_Rset(i)).strong = strong( Static_Rset(i), R_set, level ); % determine the R_set points strongly connected to the point we consider
        S(Static_Rset(i)).strong = union( S(Static_Rset(i)).strong, Static_Rset(i) ); % count i in the neighborhood
        
        if ismember( S(Static_Rset(i)).strong, R_set ) % if the strongly connected neighborhood of i is contained in R_set
            neigh_count=neigh_count + 1; % count a newly aggregated neighborhood
            hoods( neigh_count ).set = []; % initialize the new neighborhood
            % add the strongly connected neighborhood of i to the new neighborhood
            hoods( neigh_count ).set = union( hoods( neigh_count ).set, S(Static_Rset(i)).strong );
            % remove the strongly connected neighborhood of i from the R_set
            R_set = setdiff(R_set, hoods( neigh_count ).set);
        end
    % else
        % divya  insert zeros here
        % S(Static_Rset(i)).strong 
        % end divya
    end
    i=i+1;
end


% here we will attempt to add the remaining points in R_set to an existing neighborhood
Static_Rset = R_set; % make a static copy of the R_set before this step
Rsize = size(R_set); % find the initial size of the R_set (before this step)

i=1; % initialize the point consideration counter
while i <= Rsize(2) % while there are still points to be considered for neighborhood formation (in the static set)
    k=1; % initialize the neighborhood
    while k <= neigh_count % while there are still existing neighborhoods to consider
       % determine the intersect of the set of points strongly connected to the point we consider and the k-th neighborhood
       SCN = intersect( S(Static_Rset(i)).strong, hoods( k ).set );
       if size(SCN) % if that intersect is non-empty
           hoods( k ).set = union( hoods( k ).set, Static_Rset(i) ); % add i to neighborhood k
           R_set = setdiff(R_set, Static_Rset(i)); % remove i from the R_set
           k=neigh_count+1; % end the neighborhood search loop
       end
       k=k+1;
    end
   i=i+1;
end

% here we will make all the remaining points in R_set into aggregates consisting of strongly coupled sub-neighborhoods
% NOTE: this step is typically not performed since R_set is usually empty after the previous step
Static_Rset = R_set; % make a static copy of the R_set before this step
Rsize = size(R_set); % find the initial size of the R_set (before this step)

i=1; % initialize the R_set counter
while i <= Rsize(2) % while there are still points to be considered for neighborhood formation (in the static set)
    neigh_count=neigh_count + 1; % count a newly aggregated neighborhood
    hoods( neigh_count ).set = []; % initialize the new neighborhood
    % determine the intersect of the set of points strongly connected to i and the R_set
    SCR = intersect( S(Static_Rset(i)).strong, R_set );
    % add that intersected set to the new neighborhood
    hoods( neigh_count ).set = union( hoods( neigh_count ).set, SCR );
    % remove that intersected set from the R_set
    R_set = setdiff(R_set, SCR);
    i=i+1;
end

% count = neigh_count
% remaining=size(R_set)
% for counter=1:neigh_count
%     neighborhood = hoods(counter).set
% end

% ----------------------- SMOOTHING OF AGGREGATION -----------------------% 

% define the tentative prolongation P_l
Rsize = size(Original_Rset);
% P_l(i,j) = 1 if node i is in the jth neighborhood, zero otherwise
P_l = sparse(zeros(Rsize(2),neigh_count)); % Sparse added bt Derrick Cerwinsky
for i=1:Rsize(2) % for all points in the original Rset
    for j=1:neigh_count % iterate through all neighborhoods
        [is,dummy] = ismember(Original_Rset(i), hoods( j ).set); % check if the ith point is in the jth hood
        if is % if it is
            P_l(i,j) = 1; % set the prolongation operator at tindex i,j to 1
        end
    end
end

% improve the tentative prolongation via a simple smoothing step
% for this version I have hard coded a damped jacobi with parameter

% determine the "filtered" fine stiffness
AFILT = zeros(nnodes,nnodes);
% Divya----------------------------------------------------------
%  flag = 0;
% ---------------------------------------------------------------
for i=1:nnodes % for all points in the fine stiffness matrix
    for j=1:nnodes
%          % Divya--------------------------------------------------
%          if((size(S,2) == 49) && (flag == 0))
%              i
%              j
%              nnodes
%          end
%          if((size(S,2) == 20) && (flag == 0))
%              i
%              j
%              nnodes
%              flag = 1;
%          end
%          % -----------------------------------------------------
        if(i~=j) % off the diagonal of the stiffness matrix
            
            %  divya added if -else
            if(i<= size(S))
            [is,dummy] = ismember(j, S(i).strong); % check if the jth point is strongly connected to the ith point
            else
                is = 0;
            end
            if is % if it is 
                AFILT(i,j) = A(level).matrix(i,j); % copy the ij entry from the fine stiffness matrix to the filtered
            else % otherwise
                AFILT(i,j) = 0; % set it to zero
            end
        else % when you are ON the diagonal of the stiffness matrix
            sigma = 0;
            % find the sum term as described in part 2, step 11, Mandel
            for(k=1:nnodes)
                if(k~=i)
                    sigma = sigma + ( A(level).matrix(i,k) - AFILT(i,k) );
                end
            end
            AFILT(i,i) = A(level).matrix(i,i) - sigma; % add that summation as the entry in the diagonal of the filtered matrix
        end
    end
end

I = eye(Rsize(2),Rsize(2)); % identity matrix (square, dimension is number of course points)
dparam = 0.5; % damping parameter for jacobi transformation
DAMP = diag(diag(AFILT)); % get the diagonal of A
DAMP = -1 * dparam * inv(DAMP); % get the negative damped diagonal inverse of A

JACMAT = I-DAMP*AFILT; % damped jacobi transformation matrix (I - (dparam D^-1) A)

% -------------------- CALCULATE INTERPOLATION WEIGHTS --------------------% 

W(level).Iweight = JACMAT * P_l; % apply the damped jacobi transformation matrix to the initial prolongation

% -------------------------- COARSE LEVEL SETUP --------------------------% 

% generate coarse grid matrix via the Galerkin approach
W(level).Rweight = W(level).Iweight'; % restriction is the transpose of the interpolation operator in this case
A(level+1).matrix = sparse(W(level).Rweight * A(level).matrix * W(level).Iweight); % Sparse added by Derrick Cerwinsky

% I am assuming that SA coarsening on FD problems cannot go too deep...
toodeep = 0;