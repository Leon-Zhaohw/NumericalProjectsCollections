function [toodeep] = AMGeCoarsen(level, res)
% toodeep will be set to 1 if the next coarse level could not be set up, otherwise

% Divya Bansal
% Department of Computational Sciences
% University of Kentucky

% Uses subroutines written by Gundolf Haase
% as an input to AMGeCoarsen(level, residual) we take two arguments
% res is the residual at current level passed on from the smoother used.

amg_globals;

if DEBUG == 1
    disp('Starting AMGeCoarsen')
end

% if the edge and edge connectivity matrices are empty for this level
if isempty(A(level).edges ) % determine them...
    % derive node numbering for edge matrices from element connectivity
    [Ee,ea] = GetEdgeElements_Schur( A(level).elements , A(level).elconn , PAIRS);
    % Ee is edge matrix having dimensions [2,2,n_edges]
    % ea is edge connectivity matrix having dimensions [2,n_edges]
else % bring them in for use in this algorithm
    Ee = A(level).edges;
    ea = A(level).edconn;
end

% --------------------------------------- COARSENING PHASE ------------------------------------------%

% ------------------ Strength based coloring algoritm by Dr Timothy Chartier ------------------------%

omega_ = unique(A(level).elconn(:,:)); % set of all the elements in original grid

%  Determine neighborhood for nodes which are connected via an edge.
%  neigh(i).nn is the array of neighbor nodes for node i (incl. node i).
%  neigh(i).ne is the array of neighbor elements for node i.
%  neigh(i).nm is the array of neighborhood matrix for node i.

neigh = NodeNeighbourAMGe(A(level).elements , A(level).elconn);

% initialize coarse point set C_set and fine points set F_set
% C_set = [1,2,7,10,15,18,20,24];
C_set = [];
F_set = [];

nnodes = size(neigh);               % get the number of nodes in omega_

max = -1;                           % max: save value of maximum beta_(i)
max_i = 0;                          % max_i: save i for which beta_ is maximum

for j = 1 : size(omega_,1)          % for each node in the omega_
    i = omega_(j);                  % retreive node numbers from omega_ (here i is same as j)
    beta_(i)= size(neigh(i).nn,2);  % beta_(i) is degrees of freedom in the neighborhood of i.
    beta_index(i) = i;              % saves indives corresponding to beta_
end

% Sort the values ob beta_ in descending order and save their corresponding indices in another variable
[beta_, beta_index] = descendSort(beta_, beta_index);

val  = 222;
line_3_loopVar = 1;
counter = 1;
omega_transps_sze = size(omega_,1);
while line_3_loopVar == 1

    % create a set containing union of C and F
    U_set = union(C_set, F_set);

    % C_setN: set containing the coarse nodes in neighborhood only
    grtr_eq_chk = 1;

    for i = 1 : size(omega_,1)
        C_setN(i).nn = int16([]);                       % initialize
        C_setN(i).nn = intersect(C_set, neigh(i).nn);   % Ni (intersect) C
    end


    % pick i with maximal beta_(i) such that i does not belongs to C U F.
    if (counter <=size(omega_,1))
        if(ismember(beta_index(counter),U_set)==0)
            max = beta_(counter);
            max_i = beta_index(counter);
        end
    end
    counter = counter + 1;


    F_set = union(F_set, max_i);

    % order the absolute value of the entries in row i of neighborhood matrix
    % of i in descending order, including only those entries corresponding to
    % points not in Ci or C_setN

    i_row = neigh(max_i).nm(max_i,:);
    M = size(neigh(max_i).nm,2);

    q_dof = [];                 % index of q
    q_dof_val = [];             % value of q associated with the index in q_dof
    for l = 1:M,
        if (ismember(i_row(l),C_setN(max_i).nn) == 0)
            q_dof = [q_dof l];
            q_dof_val = [q_dof_val abs(i_row(l))];
        end
    end

    % q_dof is the degree of freedom associated with entry q_dof_val(v) .
    [q_dof_val, q_dof] = descendSort(q_dof_val,q_dof);

    i_rowSize = size(q_dof, 2);

    v =1;

    line_7_loopVar = 1;

    while line_7_loopVar == 1

        C_setNN = [];
        union(C_setN(max_i).nn, q_dof(v));
        C_setNN(max_i) = union(C_setN(max_i).nn, q_dof(v));

        % calculate accuracy/cost measure for both C_setN and C_setNN
        nc1 = size(C_setN(max_i).nn,2);     % number of coarse neighbor nodes of i
        nc2 = size(C_setNN(max_i),2);       % number of coarse neighbor nodes of(i + q)
        ni = size(neigh(max_i).nn,2);       % number of nodes in neighborhood of i

        gamma1 = nc1/ni;
        gamma2 = nc2/ni;


        % ----------------------------- Split A into molecules Aff, Afc, Acf and Acc  ----------------------------%

        % introduce two operators: R = [0, I] and S, where R has the dimensions of Pt and S has the dimensions of P.

        n1 = size(omega_, 1);
        n2 = n1 - nc1;
        X = zeros(nc1, n2);
        Y = eye(nc1, nc1);
        R = sparse([X, Y]);

        V = zeros(nc1, n2);
        S = sparse([Y, V]');

        % get a accumulate matrix AA calculated using all the element stiffness matrices
        % !!!!!!!!!!!!!!! insert A(level) instead
        AA = AccumulateMatrix(A(level).elements , A(level).elconn);

        Aff = sparse(S' * AA * S);
        Afc = sparse(S' * AA * R');
        Acf = sparse(R * AA * S);
        Acc = sparse(R * AA * R');

        % ---------- Calculate Interpolation operator P = [Pfc,Ic]' where Pfc = -(Afc' * inv(Aff))------------------%
        Pfc = -(Afc' * inv(Aff));
        U = zeros(nc1,(n1 - (size(Pfc,1) - size(Y,1))));
        W(level).Iweight = [Pfc, Y, U]';                % W(level).Iweight is P and Y = Ic, calculated above
        W(level).Rweight = W(level).Iweight';           % restriction is the transpose of the interpolation operator in this case

        [r,s] = size(Pfc);
        p = size(Y,1);
        u = zeros(r,n2);
        w = zeros(p,n2);

        T = n1-(size(u,1) - size(w,1));                 % to fill in Q for dimension n1*n1
        Q = sparse([u Pfc; w Y; zeros(T, n1)]);

        Basis_vector = [];


        % ---------------------------------------- Construct local Q -----------------------------------------------%
        % for constructing local Q ie q, consider its rows corresponding to each degree of freedom i belongs to F
        for i = 1 : size(F_set,2)
            if F_set(i) == max_i
                Q_for_q(max_i).nn = Q(F_set(i), :);
            end
            Q_for_q(i).nn = Q(F_set(i), :);
        end

        nQ_for_q = size(Q_for_q(1).nn,2);               % number of columns in Q_for_q
        Basis_vector = eye(nQ_for_q);

        for i = 1 : size(F_set,2)
            % ?????????????????  put check that q belongs to Zi from page 28
            if F_set(i) == max_i
                q(max_i).nn = [];                       % q is local Q
                q(max_i).nn = Basis_vector(:,max_i) * Q_for_q(max_i).nn;
            end
            q(i).nn = [];                               % q is local Q
            q(i).nn = Basis_vector(:,i) * Q_for_q(i).nn;
        end


        % ------------------- Find eigenvalue e of the stiffness matrix ------------------%
        [e,x] = eigs(neigh(max_i).nm);   % put check for e ie e does not belong to Null(Bi) page 28??????????????????????


        % ----------------------- calculate interpolation quality measure Mi --------------%
        Bi_max_i = neigh(max_i).nm(max_i,:);
        M_iq = IQMeasure(e, Bi_max_i, q(max_i).nn, Basis_vector);


        % -------------------- calculate accuracy/cost measure --------------------- %
        u1 = ACMeasure(M_iq, neigh(max_i).nm, gamma1);
        u2 = ACMeasure(M_iq, neigh(max_i).nm, gamma2);

        % code step 7 and 8 of algorithm
        if u2 >= u1
            C_set = union (C_set, q_dof(v));
            v = v+1;

            for j = 1:i_rowSize
                if(size(intersect(q_dof(v), neigh(q_dof(j)).nn),2) ~= 0)
                    beta_(j) = beta_(j) + 1;
                end
            end

        else
            line_7_loopVar = 0;
        end

    end % end of while loop

    union_temp = union(C_set,F_set);

    if size(intersect(union_temp,omega_'),2) == omega_transps_sze
        line_3_loopVar = 0;
    end

end % end of while loop



% -------------------------- COARSE LEVEL SETUP --------------------------%


% determine the number of course edges which were selected
% numcourseedges = size(A(level+1).edconn,2);
% if( numcourseedges < 3 ) % not enough course edges are selected to form an element
% NOTE: it would also be good to make sure that the 3 or more edges actually
% form an element. I propose using a simple Kruskhal's algorithm to
% determine if a circular link is present between the existing edges
%    toodeep = 1;
% else
toodeep = 0;
% end

%For now, set the next level to equal this level
A(level+1).elements = A(level).elements;
A(level+1).elconn = A(level).elconn;

% generate coarse grid matrix via the Galerkin approach
A(level+1).matrix = sparse(W(level).Rweight * A(level).matrix * W(level).Iweight);

if DEBUG == 1
    disp('Exit AMGeCoarsen')
end