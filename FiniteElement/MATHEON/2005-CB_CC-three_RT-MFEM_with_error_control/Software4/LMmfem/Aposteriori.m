
function eta_T = Aposteriori(element,coordinate,...
                 dirichlet,Neumann,u,pEval)
             
u_h = zeros(size(coordinate,1),2);
supp_area = zeros(size(coordinate,1),1);
WeightedPoint = reshape(sum(reshape(coordinate(element',:), ...
                              3,2*size(element,1))),size(element,1),2)/3;                                     
for j = 1 : size(element,1)
   supp_area(element(j,:)) = supp_area(element(j,:)) + ...
   ones(3,1)*det([1,1,1;coordinate(element(j,:),:)'])/6;
   u_h(element(j,:),:) = u_h(element(j,:),:) + ...
         det([1,1,1;coordinate(element(j,:),:)']) * ...
         ( (pEval(3*(j-1)+[1,2,3],:))' * [4 1 1;1 4 1;1 1 4]/36)';
end

temp_uh=u_h./(supp_area*ones(1,2));
u_h=tangent(coordinate,dirichlet,temp_uh);

eta_T=zeros(size(element,1),1);
for j = 1 : size(element,1)  
  eta_T(j) = sqrt(det([1,1,1;coordinate(element(j,:),:)']) * ...
              (sum( ( [4 1 1;1 4 1;1 1 4]/6 * u_h(element(j,:),:) -...
               pEval(3*(j-1)+[1,2,3],:)).^2') * ones(3,1)/6 ));           
end        

eta_T = sqrt(sum(eta_T.^2));  
 

















