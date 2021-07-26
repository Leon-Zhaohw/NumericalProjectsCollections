function [Ee,ea] = GetEdgeElements_Schur(Ae,ia,pairs)
% Compute the edge element matrices [J. Kraus] via Schur complement. 
% GetEdgeElements_Schur(Ae,ia,pairs) has input parameters 
% element matrices Ae[maxvertex,maxvertex,nelem], 
% element connectivity ia[maxvertex,nelem] 
% (local to global node numbering),
% description of edges via local node numbers in the 
% chosen class of elements via pairs[maxvertex,2].
%
% Output parameters are
% edge element matrices Ee[2,2,nedges],
% edge connectivity matrix ea[2,nedge] (local to global edge numbering).


% l2g_edges[maxedge,nelem] (local to global mapping of element edges)

amg_globals;

if DEBUG == 1
    disp('Starting GetEdgeElements_Schur')
end

% Check whether ea is really sorted.
nelem = length(Ae);
maxedge = length(pairs);        % number of edges per element
maxvertex = max(max(pairs));    % number of vertices per element
idx = int16(1:maxvertex);

%   derive node numbering for edge matrices from element connectivity ia 
%   Only valid for linear triangles

ea = int16(zeros(nelem*maxedge,2)); % array for numering of edge nodes
l2g_edges = int16([maxedge,nelem]); % local to global edge numbering per element
ik = 0;
for i=1:nelem
    for k=1:maxedge
        ik = ik+1;
        ea(ik,:) = sort([ia(pairs(k,1),i), ia(pairs(k,2),i)]);
        l2g_edges(k,i) = ik;
    end;
end;

[ea, ea_m, ea_n] = unique(ea,'rows');    %   erase double entries (=cols) in ea 
for i=1:nelem               % and renumerate also l2g_edges
    for k=1:maxedge
        l2g_edges(k,i) = ea_n(l2g_edges(k,i));
    end;
end;
ea=permute(ea,[2,1]);   % now: ea[2,nedge]
clear ea_m ea_n;

%   Compute the edge element matrices (scalar case)
%   Only valid for linear triangles (maxedge==3)
%
nedge = length(ea);         % number of edges
Ee = zeros(2,2,nedge);      % edge element matrices

%   now, one 'pairs' for all elements
%i2 = length(pairs(k,:));
i2 = length(pairs(maxedge,:));
i1 = maxvertex - i2;
rest = int16(zeros(i1,maxvertex));
for k=1:maxvertex
    rest(:,k) = setdiff(idx,pairs(k,:));
end;

for i=1:nelem
    for k=1:maxvertex
%         A11 = zeros(i1,i1); A12 = zeros(i1,i2); 
%         A21 = zeros(i2,i1); A22 = zeros(i2,i2);
        A11(:,:) = Ae(rest(:,k),rest(:,k),i);
        A12(:,:) = Ae(rest(:,k),pairs(k,:),i);
        A21(:,:) = Ae(pairs(k,:),rest(:,k),i);
        A22(:,:) = Ae(pairs(k,:),pairs(k,:),i);

        Ee(:,:,l2g_edges(k,i)) = Ee(:,:,l2g_edges(k,i)) + A22 - A21 * inv(A11) * A12;
    end;
end;


if DEBUG == 1
    disp('Exit GetEdgeElements_Schur')
end






%   each element may have a different 'pairs' set.
%
% for i=1:nelem
%     for k=1:maxvertex
%         rest = setdiff(idx,pairs(k,:));
%         i1 = length(rest);
%         i2 = length(pairs(k,:));
%         A11 = zeros(i1,i1); A12 = zeros(i1,i2); 
%         A21 = zeros(i2,i1); A22 = zeros(i2,i2);
%         A11(:,:) = Ae(rest,rest,i);
%         A12(:,:) = Ae(rest,pairs(k,:),i);
%         A21(:,:) = Ae(pairs(k,:),rest,i);
%         A22(:,:) = Ae(pairs(k,:),pairs(k,:),i);
% 
%         Ee(:,:,l2g_edges(k,i)) = Ee(:,:,l2g_edges(k,i)) + A22 - A21 * inv(A11) * A12;
%     end;
% end;
