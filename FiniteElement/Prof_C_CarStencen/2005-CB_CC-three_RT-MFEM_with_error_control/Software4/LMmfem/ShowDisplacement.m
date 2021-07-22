% ShowDisplacement.m file for the graphical representation for
% the approximate displacement $u$ in the paper
% "Three Matlab Implementations of the Lowest-Order Raviart-Thomas
% MFEM with a Posteriori Error Control" by
% C.Bahriawati and C.Carstensen

function ShowDisplacement(element,coordinate,u)
hold on
for j=1:size(element,1)
    trisurf([1 2 3],coordinate(element(j,:),1),coordinate(element(j,:),2),ones(3,1)*u(j)','facecolor','interp');
end    
view(-60,50);
 
