function [x,flag,relres,iter,resvec,alpha,r,coeff,Z,P]=mpcg(A,b,maxits,tol,restart,varargin)
% [x,flag,relres,iter,resvec,alpha,r,coeff,Z,P]=mpcg(A,b,maxits,tol,restart,M1,...)
%
% A is the matrix to solve, b the right-hand side, maxits is the maximum number of
% iterations to take, tol is the relative residual tolerance for convergence,
% restart==0 means use full MPCG and >=1 means use truncated form (and
% controls how far back to look).
%
% The supplied preconditioners M1, etc. can be matrices, in which case they
% are applied with the backslash operator, e.g. M1\r, or they can be
% function handles in which case they are simply called on r.
% They can also be a cell array of preconditioners.
%
% x is the final solution, flag indicates convergence or not, relres gives
% the relative residual of the final solution, iter gives the iteration count,
% resvec gives the residual norm history, r gives the actual residual vectors,
% coeff gives the norms of the orthogonalization coefficient matrices, Z gives
% the set of all preconditioned residuals, P gives the set of all search directions,
% and alpha gives the history of alphas.

% TODO: run a compact version if not all arguments output

flag=1;     % It ain't over until the fat lady sings
n=length(b);
if isa(varargin{1},'cell')
   k=length(varargin{1});   % number of preconditioners
   if length(varargin)~=1
      error('if you pass preconditioners in a cell array you are only allowed one cell array');
   end
else
   k=length(varargin);  % number of preconditioners
end

% first step - steepest descent
x=zeros(n,1);
r=b;
Z=multiprecondition(varargin,r);
P=Z;
alpha=pinv(P'*A*P)*(P'*r);
x(:,2)=P*alpha;
r(:,2)=b-A*x(:,2);

% then MPCG steps
for i=2:maxits
  istart=k*(i-1)+1;
  iend=k*i;
  Z(:,istart:iend)=multiprecondition(varargin,r(:,i));
  P(:,istart:iend)=Z(:,istart:iend);
  if restart
    jstart=max(1,i-restart);
  else
    jstart=1;
  end
  for j=jstart:i-1
    Pj=P(:,k*(j-1)+1:k*j);
    c=Pj'*A*Z(:,istart:iend);
    coeff(i-1,j)=norm(c,'fro');
    P(:,istart:iend)=P(:,istart:iend)-Pj*pinv(Pj'*A*Pj)*c;
  end
  Pnew=P(:,istart:iend);
  alpha(:,i)=pinv(Pnew'*A*Pnew)*(Pnew'*r(:,i));
  x(:,i+1)=x(:,i)+Pnew*alpha(:,i);
  r(:,i+1)=r(:,i)-A*Pnew*alpha(:,i);
  if norm(r(:,i+1))<tol*norm(b)
    flag=0;
    break;
  end
end

resvec=sqrt(sum(r.^2));
iter=length(resvec);
relres=resvec(end)/norm(b);

% apply multiple preconditioners to the residual r
function z=multiprecondition(pre,r)

if isa(pre{1},'cell')
   z=multiprecondition(pre{1},r);
else
   z=zeros(length(r),length(pre));
   for i=1:length(pre)
      if isa(pre{i},'function_handle')
         z(:,i)=pre{i}(r);
      else
         z(:,i)=pre{i}\r;
      end
   end
end