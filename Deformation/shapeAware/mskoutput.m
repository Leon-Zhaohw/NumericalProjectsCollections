function output = mskoutput(res)
% Used by the MOSEK compability toolkit.

if ( isfield(res,'info') )
   output.iterations = res.info.MSK_IINF_INTPNT_ITER;
else
   output.iterations = 0;
end   

output.algorithm = 'large-scale: interior-point';