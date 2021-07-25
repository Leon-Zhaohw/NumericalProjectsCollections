function [SBNengine,a0marg_Jtree,damarg_Jtree, atmarg_Jtree, Mmarg_Jtree]=ConstCrack_SBNModel(T,ev,a0_CPD,da_CPD,at_CPD,M_CPD)


%Build the Bayesian Network
NN=T*3+1;

% Define the linkage relation
dag=zeros(NN);
dag(1,3)=1;
for i=1:T
    dag(i*3,i*3+1)=1;
    dag(i*3-1,i*3)=1;
end
for i=1:T-1
    dag(i*3,i*3+3)=1;
end

% Define the nodal size
for i=1:T
   Lda(i)=length(da_CPD{i});
   M_SZ=size(M_CPD{i});
   Lat(i)=M_SZ(1);
   LM(i)=M_SZ(2);
end
%
ns=zeros(1,NN);
ns(1)=length(a0_CPD);

for i=1:T
    ns(i*3-1)=Lda(i);
    ns(i*3)=Lat(i);
    ns(i*3+1)=LM(i);
end
dnodes=1:NN;
bnet=mk_bnet(dag,ns,'discrete',dnodes);

bnet.CPD{1}=tabular_CPD(bnet,1,a0_CPD);
for i=1:T
    bnet.CPD{i*3-1}=tabular_CPD(bnet,i*3-1,da_CPD{i}); 
    bnet.CPD{i*3}=tabular_CPD(bnet,i*3,at_CPD{i});
    bnet.CPD{i*3+1}=tabular_CPD(bnet,i*3+1,M_CPD{i});
end


% Inference Engine
SBNengine = jtree_inf_engine(bnet);

[SBNengine_ev,loglik_ev]=enter_evidence(SBNengine,ev);
% Generate the posterior distribution and store it
marg_ev=marginal_nodes(SBNengine_ev,1);
a0marg_Jtree=marg_ev.T;
for t=1:T
    marg_ev=marginal_nodes(SBNengine_ev,t*3-1);
    damarg_Jtree{t}=marg_ev.T;
    marg_ev=marginal_nodes(SBNengine_ev,t*3);
    atmarg_Jtree{t}=marg_ev.T;
    marg_ev=marginal_nodes(SBNengine_ev,t*3+1);
    if length(marg_ev.T)==1
        marg_ev=marginal_nodes(SBNengine_ev,t*3+4);
        Marg_Jtree{t}=marg_ev.T;
    end
    Mmarg_Jtree{t}=marg_ev.T;
end


