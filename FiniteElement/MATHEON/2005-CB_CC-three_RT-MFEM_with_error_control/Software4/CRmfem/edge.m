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
interioredge=[];
edge2element=zeros(size(i,1),4);
for m = 1:size(element,1)
  for k = 1:3
    initial_edge = element(m,k);
    end_edge   = element(m,rem(k,3)+1);
    p = nodes2edge(initial_edge,end_edge);
   %if edge2element(p,1)~=nodes2edge(initial_edge,end_edge)
   if edge2element(p,1)==0  
        edge2element(p,:)=[initial_edge end_edge nodes2element(initial_edge,end_edge) ...
                                 nodes2element(end_edge,initial_edge)]; 
        if nodes2element(end_edge,initial_edge)~=0
           interioredge(size(interioredge,1)+1,:) = [initial_edge end_edge nodes2element(initial_edge,end_edge) ...
                                              nodes2element(end_edge,initial_edge)];                                         
        end                                      
    end
  end
end
  

