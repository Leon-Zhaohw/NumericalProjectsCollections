function [X_apx, iterations] = jacobi_smoother(level, b, u_in)

% Ryan McKenzie
% Department of Computational Sciences
% University of Kentucky

% Derrick Cerwinsky
% Changing matrixes to sparse

amg_globals;
max_iter = ABSOLUTE_MAX_ITRS;

iD = sparse(diag((diag(A(level).matrix)).^(-1)));
L = sparse(-tril(A(level).matrix,-1));
U = sparse(-triu(A(level).matrix,1));
normb = norm(b);
iter = max_iter;
x=u_in;
er = [];
i = 1;
while (i <= max_iter)
  x = x - iD*(A(level).matrix*x - b);
  if (normb)
    er = [er; norm(A(level).matrix*x - b)/normb];
  else
    er = [er; norm(x)];
  end
  if (STOP_TYPE==RESID_THRESHOLD)
      if( er(i) <= STOP_VALUE )
          break;
      end
  elseif (STOP_TYPE==RESID_REDUCE)
      if( er(i) <= (er(1) * STOP_VALUE) )
          break;
      end
  elseif (STOP_TYPE==MAX_ITRS)
      if( i >= STOP_VALUE )
          break;
      end      
  end
  i = i + 1;
end
iterations = i;
X_apx = x;