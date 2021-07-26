function mskstatus(callproc,verb,dual,rcode,res)
% Internal function used by linprog, quadprog, etc.
%
%% Copyright (c) 1998-2012 MOSEK ApS, Denmark. All rights reserved.

if ( isfield(res,'symbcon') )
    sc = res.symbcon;
else    
    [r,res2] = mosekopt('symbcon');
    sc      = res2.symbcon;
end

switch ( rcode )
case { sc.MSK_RES_ERR_INV_PROBLEM }
   disp([callproc ': Invalid problem which MOSEK cannot handle.']);
   disp([callproc ': Probably the problem is nonconvex.']);
 case { sc.MSK_RES_ERR_MISSING_LICENSE_FILE, sc.MSK_RES_ERR_LICENSE,sc.MSK_RES_ERR_LICENSE_EXPIRED,sc.MSK_RES_ERR_LICENSE_VERSION,sc.MSK_RES_ERR_SIZE_LICENSE ,sc.MSK_RES_ERR_PROB_LICENSE,sc.MSK_RES_ERR_FILE_LICENSE,sc.MSK_RES_ERR_MISSING_LICENSE_FILE,sc.MSK_RES_ERR_SIZE_LICENSE_CON,sc.MSK_RES_ERR_SIZE_LICENSE_VAR,sc.MSK_RES_ERR_SIZE_LICENSE_INTVAR,sc.MSK_RES_ERR_OPTIMIZER_LICENSE,sc.MSK_RES_ERR_FLEXLM,sc.MSK_RES_ERR_LICENSE_SERVER,sc.MSK_RES_ERR_LICENSE_MAX,sc.MSK_RES_ERR_LICENSE_FEATURE,sc.MSK_RES_ERR_LICENSE_CANNOT_CONNECT,sc.MSK_RES_ERR_LICENSE_INVALID_HOSTID,sc.MSK_RES_ERR_LICENSE_SERVER_VERSION,sc.MSK_RES_ERR_OPEN_DL ,sc.MSK_RES_ERR_OLDER_DLL,sc.MSK_RES_ERR_NEWER_DLL}
  if ( isfield(res,'rmsg') )
    msg = sprintf('%s : %s',callproc,res.rmsg);
    disp(msg);
  end
end
  



if verb > 0
   if ( isfield(res,'sol') )
      if ( ( ~dual & res.sol.itr.prosta==sc.MSK_PRO_STA_PRIM_INFEAS ) | ...
           ( dual & res.sol.itr.prosta==sc.MSK_PRO_STA_PRIM_AND_DUAL_INFEAS ) )
         disp('The primal problem is infeasible.');
      end
      if ( ( ~dual & res.sol.itr.prosta==sc.MSK_PRO_STA_DUAL_INFEAS ) | ...
           ( dual & res.sol.itr.prosta==sc.MSK_PRO_STA_PRIM_INFEAS ) )
         disp('The dual problem is infeasible.');
      end
   end   
   
   if ( mskeflag(rcode,res)==1 )
      disp('Optimization terminated successfully.');   
   end
   
   if ( mskeflag(rcode,res)==0 )
     disp('Maximum number of iterations exceeded.')
   end
end
