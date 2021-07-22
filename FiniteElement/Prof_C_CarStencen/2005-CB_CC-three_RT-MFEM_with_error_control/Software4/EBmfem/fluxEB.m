% FluxEB.m file for the approximate flux for EBmfem in the paper
% "Three Matlab Implementations of the Lowest-Order Raviart-Thomas
% MFEM with a Posteriori Error Control" by
% C.Bahriawati and C.Carstensen
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The rows of (1) and (2) form systems of
%equations with matrix $A$ and the right-hand side $b=(b_D,b_f)^T$ 
%and $b=(b_D,b_f,0,b_g)^T$ for the edge
%basis function and the Lagrange multiplier technique, respectively. It is
%solved by the binary Matlab operator $\backslash$ and reads
%  x = A \ b;
%The solution vectors are $x=(x_{\psi},x_{u})$ and
%$x=(x_{\psi},x_{u},x_{\lambda_M},x_{\lambda_N})$, where $p$
%and $u$ occupy the first component and the second component for the
%edge-basis function and the Lagrange multiplier technique, respectively.
%The displacement and flux $u,p$ are extracted by
%\texttt{x(noedges+1:noedges\\+size(element,1))} and \texttt{x(1:noedges)} for
%the edge-basis function, and \texttt{x(3*size(ele\\ment,1)+1:4*size(element,1))}
%and \texttt{x(1:3*size(element,1))} for the Lagrange multiplier technique.
%Then, for each $T\in\T$, the flux $p$ is computed by
%\begin{equation}\label{eqD.1}
%  p_{T_{EB}} = \sum_{j=1}^3\psi_{E_j}({\bf x}-P_j)\,p_{EB_j},\qquad  p_{T_{LM}} = \sum_{j=1}^3\psi_{j}\,p_{LM_j},
%\end{equation}
%where $x$ belongs to the triangle $T$ with vertices $P_1$, $P_2$, and  $P_3$.
    
function p=fluxEB(element,coordinate,u,noedges,nodes2edge,edge2element)
p=zeros(3*size(element,1),2);
for j=1:size(element,1)
  signum=ones(1,3);
  signum(find(j==edge2element(diag(nodes2edge(element(j,[2 3 1]),...
  element(j,[3 1 2]))),4)))=-1;
  c=coordinate(element(j,[2 3 1]),:)-coordinate(element(j,[3 1 2]),:);
  n=[norm(c(1,:)),norm(c(2,:)),norm(c(3,:))];
  coord=coordinate(element(j,:),:)';
  N=coord(:)*ones(1,3)-repmat(coord,3,1);
  pc=reshape(N*diag(signum)*diag(n)*u(diag(nodes2edge(element(j,[2 3 1]),...
             element(j,[3 1 2]))))/det([1 1 1;coordinate(element(j,:),:)']),2,3);
  p(3*(j-1)+[1,2,3],:)=[pc(1,:)',pc(2,:)'];
end


