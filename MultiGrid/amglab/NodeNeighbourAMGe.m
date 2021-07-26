function neigh = NodeNeighbourAMGe(Ae, ia)

% Divya Bansal
% Department of Computational Sciences
% University of Kentucky

amg_globals;

% NodeNeighbourAMG(Ae,ia) will calculate the neighborhood matrix. It has
% input parameters:
% element matrices Ae[maxvertex,maxvertex,nelem], 
% element connectivity ia[maxvertex,nelem] 


% canonical basis vectors with dimension 1x3 because element stiffness
% matrix is 3x3

% Edit by Derrick Cerwinsky

%if issparse(ia) == false    % Test if ia is sparse, make sparse if not.
%    ia = sparse(ia);
%end

%if issparse(Ae) == false    % Test if Ae is sparse, make sparse if not.
%    Ae = sparse(Ae);
%end

% End edit




nnode = max(max(ia(:,:)));      % number of nodes
nelem = size(ia,2);             % number of elements 

BASIS_VEC = eye(nnode); % canonical basis vectors for all nodes in the grid
n = size(BASIS_VEC,1);

for i=1:nnode
    neigh(i).nn = int16([]);
    neigh(i).ne = int16([]);
    neigh(i).nm = int16([]);
end;

% -------------------------------------------------------------------------
% calculate the degree of freedom for each element
for i=1:nelem; % for all the elements in the grid
    
    % code to accumulate an element stiffness matrix into a matrix of
    % dimension equal to the main matrix A
    KK = sparse(zeros(nnode));
    for k=1:size(Ae(:,:,i),1)
        ik = ia(k,i);
        for l=1:size(Ae(:,:,i),2)
            il = ia(l,i);
            KK(ik,il) = KK(ik,il) + Ae(k,l,i);
        end;
    end;
    
    k = 1;
    for j=1:n, % for all nodes in a grid
        temp = BASIS_VEC(j,:) * KK * BASIS_VEC(j,:)';
        if temp ~= 0
            neigh(i).dof(k) = j;
            k = k + 1; % if temp is satisfied then add temp to dof array 
                       % and increment array count
        end
    end
    % neigh(i).dof
end

% -------------------------------------------------------------------------
% set of elements in the neighborhood of node i
for i=1:nnode, % for all nodes in the grid
    m =1;
    for j=1:nelem, % for all elements in the grid
        
        KK = sparse(zeros(nnode));
        for k=1:size(Ae(:,:,j),1)
            jk = ia(k,j);
            for l=1:size(Ae(:,:,j),2)
                jl = ia(l,j);
                KK(jk,jl) = KK(jk,jl) + Ae(k,l,j);
            end;
        end;
        
        temp = BASIS_VEC(i,:) * KK * BASIS_VEC(i,:)';
        if temp ~= 0
            neigh(i).ne(m) = j;
            m = m + 1;
        end    
    end
    % neigh(i).ne
end

% -------------------------------------------------------------------------

% set of points in the neighborhood of node i
for i=1:nnode
    i_elem = neigh(i).ne; % retrieve elements attached to node i  
    indx = size(i_elem,2); % get the size if vectot i_elem
    for j = 1:indx 
        % add the dof element-wise which are in neighbourhood of node i
        neigh(i).nn = union(neigh(i).nn, neigh(i_elem(j)).dof); 
    end
    % neigh(i).nn
end

% -------------------------------------------------------------------------

% local matrix neigh.nm on node i
X = int16([]);
Y = int16([]);
for j =1:nnode,
    i_elem = neigh(j).ne; % retrieve elements attached to node i  
    indx = size(i_elem,2); % get the size if vector i_elem
    KK = sparse(zeros(nnode));
    for i = 1:indx
        for k=1:size(Ae(:,:,i_elem(i)),1)
            ik = ia(k,i_elem(i));
            for l=1:size(Ae(:,:,i_elem(i)),2)
                il = ia(l,i_elem(i));
                KK(ik,il) = KK(ik,il) + Ae(k,l,i_elem(i));
            end;
        end;
    end;
    neigh(j).nm = sparse(KK);
%    neigh(j).nm 
end;