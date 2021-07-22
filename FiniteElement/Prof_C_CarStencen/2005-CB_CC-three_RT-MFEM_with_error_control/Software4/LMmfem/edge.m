% File edge.m

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
    p = nodes2edge(element(m,k),element(m,rem(k,3)+1)); 
    if edge2element(p,1)==0  
      edge2element(p,:)=[element(m,k) element(m,rem(k,3)+1)  nodes2element(element(m,k),element(m,rem(k,3)+1)) ...
                          nodes2element(element(m,rem(k,3)+1),element(m,k))];            
               
    end
  end
end
  
interioredge=edge2element(find(edge2element(:,4)),:);
exterioredge=edge2element(find(edge2element(:,4)==0),[1,2,3]);

