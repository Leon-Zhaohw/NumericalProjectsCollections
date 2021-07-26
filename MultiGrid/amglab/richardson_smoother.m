function [X_apx, iterations] = richardson_smoother(level, b, u_in)

% Divya Bansal
% Department of Computer Science
% University of Kentucky

% Derrick Cerwinsky
% Added sparse

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
% eigenvalues = 1:49
 

 eigenvalues = eig(A(level).matrix);
 
 i_count =1;
len = size(eigenvalues);
while (i_count<= len)
    if (eigenvalues(i_count)<0)
        eigenvalues(i_count)= eigenvalues(i_count) * -1;
    end
    i_count = i_count + 1;
end

spec_rad = max(eigenvalues);

% spec_rad = spectral_radius(level); 

while (i <= max_iter)
  x = x + spec_rad*(b - A(level).matrix * x);
  
  % message = 'Richardson'
  % x'
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