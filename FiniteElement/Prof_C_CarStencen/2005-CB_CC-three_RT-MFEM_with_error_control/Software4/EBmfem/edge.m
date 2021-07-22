%function [nodes2element,inneredge]=edge(element,coordinate);
function [nodes2element,nodes2edge,noedges,edge2element,interioredge]=edge(element,coordinate);
% to enumerate number of edges
nodes2element=sparse(size(coordinate,1),size(coordinate,1));
for j=1:size(element,1)
    nodes2element(element(j,:),element(j,[2 3 1]))=...
      nodes2element(element(j,:),element(j,[2 3 1]))+j*eye(3,3);
end

B=nodes2element+nodes2element';
[i,j]=find(triu(B));
nodes2edge=sparse(i,j,1:size(i,1),size(coordinate,1),size(coordinate,1));
nodes2edge=nodes2edge+nodes2edge';
noedges=size(i,1);

% to generate element of edge
edge2element=zeros(size(i,1),4);
for m = 1:size(element,1)
  for k = 1:3
    %initial_edge = element(m,k);   
    %end_edge = element(m,rem(k,3)+1);
    p = nodes2edge(element(m,k),element(m,rem(k,3)+1)); % update on 13/2/3
    %p = nodes2edge(initial_edge,end_edge);
    if edge2element(p,1)==0  
       %edge2element(p,:)=[initial_edge end_edge nodes2element(initial_edge,end_edge) ...
        %                  nodes2element(end_edge,initial_edge)]; 
      edge2element(p,:)=[element(m,k) element(m,rem(k,3)+1)  nodes2element(element(m,k),element(m,rem(k,3)+1)) ...
                          nodes2element(element(m,rem(k,3)+1),element(m,k))];            
               
    end
  end
end
  
% To produce the interior edges 
%[l,m]=find(edge2element(:,4));
interioredge=edge2element(find(edge2element(:,4)),:);


%[r,s]=find(edge2element(:,4)==0);
exterioredge=edge2element(find(edge2element(:,4)==0),[1,2,3]);

 


% ===============================

%edge2element1=zeros(noedges ,4);
%[I,J]=find(nodes2element);Z=[I,J];
%[K,P,Q]=unique(sort([I,J]')','rows')
%[L,M]=setdiff(Z,Z(P,:),'rows');
%edge2element1(:,1:3)=[Z(P,:) ,diag(nodes2element(Z(P,1), Z(P,2) )) ];
%edge2element1(Q(M),4)=diag(nodes2element( L(:,1),L(:,2)) );
%edge2element1

