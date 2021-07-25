function [SBNengine,damarg_Jtree, atmarg_Jtree, Mmarg_Jtree]=ConstCrack_JtreeModel(T,ev,da_CPD,a1_CPD,at_CPD,M_CPD)

Lda=length(da_CPD);
Lat=length(M_CPD(:,1));
LM=length(M_CPD(1,:));

%Build the Bayesian Network
NN=T*3;
hnodes1=1:3:3*T;
hnodes2=2:3:3*T;
onodes=3:3:3*T;

eclass=ones(1,NN);

H1class=1;
H2class=2;
H3class=4;
O1class=3;

eclass(hnodes1)=H1class;
eclass(hnodes2(1))=H2class;
eclass(hnodes2(2:end))=H3class;
eclass(onodes(1:end))=O1class;
% Define the linkage relation
dag=zeros(NN);
for i=1:T
    dag(i*3-1,i*3)=1;
    dag(i*3-2,i*3-1)=1;
end
for i=1:T-1
    dag(i*3-1,i*3+2)=1;
end
% Define the nodal size
ns=zeros(NN,1);
for i=1:3:NN
    ns(i)=Lda;
    ns(i+1)=Lat;
    ns(i+2)=LM;
end

dnodes=1:NN;
bnet=mk_bnet(dag,ns,'discrete',dnodes,'equiv_class',eclass);

bnet.CPD{H1class}=tabular_CPD(bnet,hnodes1(1),da_CPD);

bnet.CPD{H2class}=tabular_CPD(bnet,hnodes2(1),a1_CPD);

bnet.CPD{H3class}=tabular_CPD(bnet,hnodes2(2),at_CPD);

bnet.CPD{O1class}=tabular_CPD(bnet,onodes(1),M_CPD);
% Inference Engine
SBNengine = jtree_inf_engine(bnet);

[SBNengine_ev,loglik_ev]=enter_evidence(SBNengine,ev);
% Generate the posterior distribution and store it
Mmarg_Jtree=zeros(LM,T);
atmarg_Jtree=zeros(Lat,T);
damarg_Jtree=zeros(Lda,T);
for t=1:T
    marg_ev=marginal_nodes(SBNengine_ev,t*3-2);
    damarg_Jtree(:,t)=marg_ev.T;
    marg_ev=marginal_nodes(SBNengine_ev,t*3-1);
    atmarg_Jtree(:,t)=marg_ev.T;
    marg_ev=marginal_nodes(SBNengine_ev,t*3);
    if length(marg_ev.T)==1
        marg_ev.T=zeros(LM,1);
    end
    Mmarg_Jtree(:,t)=marg_ev.T;
end


