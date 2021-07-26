function [x, iter]  = SOR_smoother(level, b, x )

%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
% [x, error, iter, flag]  = sor(A, x, b, w, max_it, tol)
%
% sor.m solves the linear system Ax=b using the 
% Successive Over-Relaxation Method (Gauss-Seidel method when omega = 1 ).
%
% input   A        REAL matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         w        REAL relaxation scalar
%         max_it   INTEGER maximum number of iterations
%         tol      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it

%This section is added for compatibility with AMGLab
  amg_globals;

if DEBUG == 1
    disp('Starting SOR_smoother')
end
  
  max_it = ABSOLUTE_MAX_ITRS;
%  sparse As;
  As = sparse(A(level).matrix);
  
   w = OMEGA;

  iter = 0;


  error_1 = norm( b - As*x );
 
  if(STOP_TYPE == RESID_THRESHOLD)
    if (error_1 <= STOP_VALUE)
        return;
    end
  end
  
  
%  sparse N;
%  sparse M;
  
  % matrix splitting
  
   c = w * b;
   M =  sparse(w * tril( As, -1 ) + diag(diag( As )));
   N = sparse(-w * triu( As,  1 ) + ( 1.0 - w ) * diag(diag( As )));


    for iter = 1:max_it                         % begin iteration

        x_1 = x;
        x   = M \ ( N*x + c );                   % update approximation
        
     %   error = norm( x - x_1 ) / norm( x );     % compute error
    error = norm(b - As*x);
   
   % check convergence
      switch STOP_TYPE
          case RESID_THRESHOLD
              if (error <= STOP_VALUE)
                  break;
              end
          case RESID_REDUCE
              if( error <= error_1 * STOP_VALUE )
                  break;
              end
      end
     
    end
  
if DEBUG == 1
    disp('Exit SOR_smoother')
end
 
