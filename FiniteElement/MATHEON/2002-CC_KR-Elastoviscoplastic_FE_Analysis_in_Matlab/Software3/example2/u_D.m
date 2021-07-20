function [W,M]=u_D(x,t)
% values of the function on the dirichlet boundary
M=zeros(2*size(x,1),2);
W=zeros(2*size(x,1),1);
M(1:2:2*size(x,1)-1,1)=ones(size(x,1),1);
M(2:2:2*size(x,1),2)  =ones(size(x,1),1);
value=zeros(size(x,1),2);
W(1:2:2*size(x,1)-1,1)=value(:,1);
W(2:2:2*size(x,1),1)  =value(:,2);
