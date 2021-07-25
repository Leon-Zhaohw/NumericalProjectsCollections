function Beta=ReliabIndex(atmarg,at_R,atc)

Iatc=find(at_R==atc);
FailProb=sum(atmarg(Iatc:end));
Beta=-norminv(FailProb);

