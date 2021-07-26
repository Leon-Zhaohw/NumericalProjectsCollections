function [param,mosek_exists] = default_mosek_param()
  % DEFAULT_MOSEK_PARAM
  % 
  % Outputs:
  %   param stuct containing some nice default mosek params
  %   mosek_exists  whether mosek exists
  %

  % Tolerance parameter
  % >1e0 NONSOLUTION
  % 1e-1 artifacts in deformation
  % 1e-3 artifacts in isolines
  % 1e-4 seems safe for good looking deformations
  % 1e-8 MOSEK DEFAULT SOLUTION
  % 1e-14 smallest allowed value
  if(exist('mosekopt','file'))
    param.MSK_DPAR_INTPNT_TOL_REL_GAP = 1e-14;
    %param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 0;
    %param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 0;
    %param.MSK_DPAR_INTPNT_CO_TOL_MU_RED = 0;
    %param.MSK_DPAR_INTPNT_CO_TOL_INFEAS = 0;
    %param.MSK_IPAR_INTPNT_ORDER_METHOD = '';
    % always use one core and always leave one core
    param.MSK_IPAR_INTPNT_NUM_THREADS = max(feature('numCores')-1,1);
    if(isunix)
      % Get the real number of cores
      [r,c] = system('sysctl hw.ncpu | awk ''{print $2}''');
      if r==0
        c = str2double(c);
        if ~isnan(c)
          param.MSK_IPAR_INTPNT_NUM_THREADS = max(c-1,1);
        end
      end
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
end
