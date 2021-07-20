function [W,M]=u_D(lcoordinates,t)
M=zeros(2*size(lcoordinates,1),2);
W=zeros(2*size(lcoordinates,1),1);
% conditions on the x-axis
tmp = find(lcoordinates(:,1)>0 & lcoordinates(:,2)==0);
M(2*tmp-1,1:2) = ones(size(tmp,1),1)*[0 1];
W(2*tmp-1,1) = zeros(size(tmp,1),1);
% conditions on the y-axis
tmp = find(lcoordinates(:,2)>0 & lcoordinates(:,1)==0); 
M(2*tmp-1,1:2) = ones(size(tmp,1),1)*[1 0];
W(2*tmp-1,1) = zeros(size(tmp,1),1);
