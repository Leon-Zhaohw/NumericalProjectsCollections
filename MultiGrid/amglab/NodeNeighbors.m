function neigh = NodeNeighbors(ia)
% Determine the neighborhood of nodes.
% element connectivity ia[loc,nelem] for local to global node numbering.

nnode = max(max(ia(:,:)));      % number of nodes
nelem = length(ia);             % number of edges ( =size(ia,2) )

for i=1:nnode
    neigh(i).nn = int16([]);
end;
for i=1:nelem
    for k=1:size(ia(:,i),1)
        ik = ia(k,i);
        for l=1:size(ia(:,i),1)
            il = ia(l,i);
            neigh(ik).nn = [neigh(ik).nn il];
        end;
    end;
end;
for i=1:nnode
    neigh(i).nn = unique(neigh(i).nn);
%     neigh(i)
end;

