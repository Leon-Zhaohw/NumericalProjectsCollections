function KK = AccumulateMatrix(Ae,ia)
% Accumulates the f.e. matrices to a global f.e. matrix
% INPUT:  
% element matrices Ae[loc,loc,nelem], 
% element connectivity ia[loc,nelem] for local to global node numbering.
% The node numbering starts with 1 (Fortran style).
%

amg_globals;

if DEBUG == 1
    disp('Starting AccumulateMatrix')
end

nnode = max(max(ia));           %   number of nodes
KK = zeros(nnode,nnode);
for i=1:size(Ae,3)
    for k=1:size(Ae(:,:,i),1)
        ik = ia(k,i);
        for l=1:size(Ae(:,:,i),2)
            il = ia(l,i);
            KK(ik,il) = KK(ik,il) + Ae(k,l,i);
        end;
    end;
end;


if DEBUG == 1
    disp('Exit AccumulateMatrix')
end

