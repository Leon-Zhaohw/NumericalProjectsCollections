
function pEval=fluxLMEval(element,coordinate,u)

% ===========================================
% to compute p   by CB
% ===========================================

weightedpoint = reshape(sum(reshape(coordinate(element',:),3,2*size(element,1))) ...
                 ,size(element,1),2)/3;
P = reshape(u(1:3*size(element,1)),3,size(element,1));
pxEval=[4 1 1;1 4 1;1 1 4]/6 * reshape(coordinate(element',1),3,size(element,1));
pyEval=[4 1 1;1 4 1;1 1 4]/6 * reshape(coordinate(element',2),3,size(element,1));
p_x = ones(3,1)*P(1,:)+(ones(3,1)*P(3,:)) .* (pxEval - ones(3,1)*weightedpoint(:,1)');
p_y = ones(3,1)*P(2,:)+(ones(3,1)*P(3,:)) .* (pyEval - ones(3,1)*weightedpoint(:,2)');  
pEval=reshape([p_x,p_y],3*size(element,1),2);

