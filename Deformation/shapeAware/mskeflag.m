function [exitflag] = mskeflag(rcode,res)
% Used by the MOSEK compability toolbox.
%
%% Copyright (c) 1998-2012 MOSEK ApS, Denmark. All rights reserved.

if rcode==0 && isfield(res,'sol') && res.sol.itr.solsta==res.symbcon.MSK_SOL_STA_OPTIMAL; 
  exitflag = 1;
else
  if ( isfield(res,'symbcon') && rcode==res.symbcon.MSK_RES_TRM_MAX_ITERATIONS )
     exitflag = 0;
  else
     exitflag = -1;
  end
end  
