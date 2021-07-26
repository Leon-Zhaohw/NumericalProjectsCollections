function neigh = NodeNeighborsEdges(ia)
% Determine the neighborhood of nodes.
% element connectivity ia[loc,nelem] for local to global node numbering.

amg_globals;

if DEBUG == 1
    disp('Starting NodeNeighborsEdges')
end

nnode = max(max(ia(:,:)));      % number of nodes
nelem = size(ia,2); % length(ia);             % number of edges ( =size(ia,2) )

for i=1:nnode
    neigh(i).nn = int16([]);
    neigh(i).ne = int16([]);
end;
for i=1:nelem
    for k=1:size(ia(:,i),1)
        ik = ia(k,i);
        neigh(ik).ne = [neigh(ik).ne i];
        for l=1:size(ia(:,i),1)
            il = ia(l,i);
            neigh(ik).nn = [neigh(ik).nn il];
            neigh(il).ne = [neigh(il).ne i];
        end;
    end;
end;
for i=1:nnode
    neigh(i).nn = unique(neigh(i).nn);
    neigh(i).ne = unique(neigh(i).ne);
%     neigh(i)
end;


if DEBUG == 1
    disp('Exit NodeNeighborsEdges')
end