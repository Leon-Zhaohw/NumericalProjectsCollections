function W = biharmonic_bounded(V,F,b,bc,type,pou,low,up)
  % BIHARMONIC_BOUNDED Compute biharmonic bounded coordinates, using quadratic
  % optimizer
  %
  % W = biharmonic_bounded(V,F,b,bc,type,pou)
  % W = biharmonic_bounded(V,F,b,bc,type,pou,low,up)
  %
  % Inputs:
  %  V  list of vertex positions
  %  F  list of face indices, for 3D F is #F by 4, for 2D F is #F by 3
  %  b  list of boundary vertices
  %  bc list of boundary conditions, size(boundary) by # handles matrix of
  %    boundary conditions where bc(:,i) are the boundary conditions for the 
  %    ith handle
  %  Optional:
  %    type  type of optimizer to use {best available}:
  %      'quad'
  %      'least-squares'
  %      'conic'
  %    pou  true or false, enforce partition of unity explicitly {false}
  %    low  lower bound {0}
  %    up  upper bound {1}
  %  
  %
  % Outputs:
  %  W  weights, # vertices by # handles matrix of weights
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: boundary_conditions
  %

  % set default for enforcing partition of unity constraint
  if ~exist('pou','var') || isempty(pou)
    pou = false;
  end

  % number of vertices
  n = size(V,1);
  % number of handles
  m = size(bc,2);

  % Build discrete laplacian and mass matrices used by all handles' solves
  if(size(F,2)==4)
    fprintf('Solving over volume...\n');
    L = cotmatrix3(V,F);
    M = massmatrix3(V,F,'barycentric');
  else
    L = cotmatrix(V,F);
    M = massmatrix(V,F,'voronoi');
  end

  % default bounds
  if ~exist('low','var') || isempty(low)
    low = 0;
  end
  if ~exist('up','var') || isempty(up)
    up = 1;
  end

  % set default optimization method
  if ~exist('type','var') || isempty(type)
    if(exist('mosekopt'))
      % if mosek is available then conic is fastest
      type = 'conic';
    else
      % if we only have matlab then quadratic is default
      type = 'quad';
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SET UP SOLVER
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % check for mosek and set its parameters
  param = [];
  % Tolerance parameter
  % >1e0 NONSOLUTION
  % 1e-1 artifacts in deformation
  % 1e-3 artifacts in isolines
  % 1e-4 seems safe for good looking deformations
  % 1e-8 MOSEK DEFAULT SOLUTION
  % 1e-14 smallest allowed value
  if(exist('mosekopt','file'))
    if(strcmp(type,'quad') || strcmp(type,'least-squares') )
      param.MSK_DPAR_INTPNT_TOL_REL_GAP = 1e-10;
    elseif(strcmp(type,'conic') )
      param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-10;
    end
    mosek_exists = true;
  else 
    mosek_exists = false;
    if(verLessThan('matlab','7.12'))
      % old matlab does not solve quadprog with sparse matrices: SLOW
      % solution: dowloand MOSEK or upgrade to 2011a or greater
      warning([ ...
        'You are using an old version of MATLAB that does not support ' ...
        'solving large, sparse quadratic programming problems. The ' ...
        'optimization will be VERY SLOW and the results will be ' ...
        'INACCURATE. Please install Mosek or upgrade to MATLAB version >= ' ...
        '2011a.']);
    else
      % Tell matlab to use interior point solver, and set tolerance
      % 1e-8 MATLAB DEFAULT SOLUTION (very low accuracy)
      % 1e-10 (low accuracy)
      % 1e-12 (medium-low accuracy)
      % 1e-14 (medium accuracy)
      % 1e-16 (high accuracy)
      param = optimset( ...
        'TolFun',1e-16, ...
        'Algorithm','interior-point-convex', ...
        ... % 'Algorithm','active-set', ...
        'MaxIter', 1000, ...
        'Display','off');
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SET UP PROBLEM AND SOLVE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(pou)
    % Enforce partition of unity as explicity constraints: solve for weights
    % of all handles simultaneously
    if(strcmp(type,'quad'))
      % biharmonic system matrix
      Qi = L*(M\L);
      Q = sparse(m*n,m*n);
      % Q is sparse matrix with Qi along diagonal
      for ii = 1:m
        d = (ii - 1)*n + 1;
        Q(d:(d + n-1), d:(d + n-1)) = Qi;
      end
      % linear constraints: partition of unity constraints and boundary
      % conditions
      PA = repmat(speye(n,n),1,m);
      Pb = ones(n,1);
      % boundary conditions
      BCAi = speye(n,n);
      BCAi = BCAi(b,:);
      BCA = sparse(m*size(BCAi,1),m*size(BCAi,2));
      % BCA is sparse matrix with BCAi along diagonal
      for ii = 1:m
        di = (ii - 1)*size(BCAi,1) + 1;
        dj = (ii - 1)*size(BCAi,2) + 1;
        BCA(di:(di + size(BCAi,1)-1), dj:(dj + size(BCAi,2)-1)) = BCAi;
      end
      BCb = bc(:);
      % set bounds
      ux = up.*ones(m*n,1);
      lx = low.*ones(m*n,1);
      if(mosek_exists)
        fprintf('Quadratic optimization using mosek...\n');
      else
        fprintf('Quadratic optimization using matlab...\n');
      end
      fprintf( [ ...
        '  minimize:     x''LM\\Lx\n' ...
        'subject to: %g <= x <= %g, ∑_i xi = 1\n'], ...
        low,up);
      tic;
      W = quadprog(Q,zeros(n*m,1),[],[],[PA;BCA],[Pb;BCb],lx,ux,[],param);
      toc
      W = reshape(W,n,m);
    else
      error( [ ...
        'Enforcing partition of unity only support in conjunction with ' ...
        'type=''quad''']);
    end
  else
    % Drop partition of unity constraints, solve for weights of each handle
    % independently then normalize to enforce partition of unity
    if(strcmp(type,'quad'))
      % build quadratic coefficient matrix (bilaplacian operator)
      Q = L*(M\L);
      % set bounds
      ux = up.*ones(n,1);
      lx = low.*ones(n,1);
    elseif(strcmp(type,'least-squares'))
      % solve same problem but as least-squares problem see mosek documention
      % for details
      I = speye(n);
      Z = sparse(n,n);
      Q = [Z,Z;Z,I];
      F = sqrt(M)\L;
      c = zeros(n,1);
      B = [F,-I];
      ux = [up.*ones(n,1) ;  Inf*ones(n,1)];
      lx = [low.*ones(n,1); -Inf*ones(n,1)];
    elseif(strcmp(type,'conic'))
      % solve same problem but as conic problem see mosek documention for
      % details
      F = sqrt(M)\L;
      prob.c = [zeros(2*n,1); 1];
      I = speye(n);
      prob.a = [F,-I,zeros(n,1)];
      prob.blc = zeros(n,1);
      prob.buc = zeros(n,1);
      prob.bux = [ up.*ones(n,1);  Inf*ones(n,1);  Inf];
      prob.blx = [ low.*ones(n,1); -Inf*ones(n,1); -Inf];
      prob.cones = cell(1,1);
      prob.cones{1}.type = 'MSK_CT_QUAD';
      t_index = 2*n +1;
      z_indices = (n+1):(2*n);
      prob.cones{1}.sub = [t_index z_indices];
    else
      error('Bad type');
    end

    % number of handles
    m = size(bc,2);
    % allocate space for weights
    W = zeros(n,m);
    tic;
    % loop over handles
    for i = 1:m
      if(strcmp(type,'quad'))
        % enforce boundary conditions via lower and upper bounds
        %lx(b) = bc(:,i);
        %ux(b) = bc(:,i);
        Aeq = speye(n,n);
        Aeq = Aeq(b,:);
        if(mosek_exists)
          fprintf('Quadratic optimization using mosek...\n');
        else
          fprintf('Quadratic optimization using matlab...\n');
        end
        fprintf( [ ...
          '  minimize:     x''LM\\Lx\n' ...
          'subject to: %g <= x <= %g\n' ], ...
          low,up);
        % if mosek is not available, then matlab will complain that sparse
        % matrices are not yet supported...
        [x,fval,err] = quadprog(Q,zeros(n,1),[],[],Aeq,bc(:,i),lx,ux,[],param);
        if(err ~= 1)
          fprintf([...
            '----------------------------------------------------------\n' ...
            'ERROR ('  num2str(err) ',' num2str(fval) '):' ...
            ' solution may be inaccurate...\n' ...
            '----------------------------------------------------------\n' ...
            ]);
        end
      elseif(strcmp(type,'least-squares'))
        % enforce boundary conditions via lower and upper bounds
        lx(b) = bc(:,i);
        ux(b) = bc(:,i);
        fprintf('Quadratic optimization using mosek...\n');
        fprintf([ ...
          '  minimize:       z''z\n' ...
          '  subject to: M\\Lx - z = 0\n' ...
          '  and          %g <= x <= %g\n'], ...
          low,up);
        x = quadprog(Q,zeros(2*n,1),[],[],B,c,lx,ux,[],param);
      elseif(strcmp(type,'conic'))
        prob.bux(b) = bc(:,i);
        prob.blx(b) = bc(:,i);
        fprintf('Conic optimization using mosek...\n');
        fprintf([ ...
          '  minimize:         t\n' ...
          '  subject to: M\\Lx - z = 0,\n' ...
          '             t >= sqrt(z''z),\n' ...
          '               %f <= x <= %f\n'], ...
          low,up);
        [r,res]=mosekopt('minimize echo(0)',prob,param);
        % check for mosek error
        if(r == 4006)
          warning(['MOSEK ERROR. rcode: ' ...
            num2str(res.rcode) ' ' ...
            res.rcodestr ' ' ...
            res.rmsg ...
            'The solution is probably OK, but ' ...
            'to make this error go away, increase: ' ...
            'MSK_DPAR_INTPNT_CO_TOL_REL_GAP' ...
            n]);
        elseif(r ~= 0)
          error(['FATAL MOSEK ERROR. rcode: ' ...
            num2str(res.rcode) ' ' ...
            res.rcodestr ' ' ...
            res.rmsg]);
        end
        % extract solution from result
        x = res.sol.itr.xx;
      end
      % set weights to solution in weight matrix
      W(:,i) = x(1:n);
      fprintf('Lap time: %gs\n',toc);
    end
    t = toc;
    fprintf('Total elapsed time: %gs\n',t);
    fprintf('Average time per handle: %gs\n',t/m);
  end

end
