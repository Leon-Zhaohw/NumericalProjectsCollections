function [atmarg_JtreeFine Beta_JtreeFine,at_RFine,Ytc_CPD,at_CPD,M_CPD]=StochasticCrack_Fine(alphaY_mean,alphaY_std,Y_mean,Y_std,a0_mean,M_sigma,DeltaN,DeltaS,m,C,atc,atLimit,ev,T)

LalphaY=100;
LY=100;
Lat=100;

alphaY_R=[0 10.^(1:(6-1)/(LalphaY-1):6)];
Y_R=[0.0:15/(LY-1):15 inf];
at_R=[0 exp(log(0.01):(log(atc)-log(0.01))/(Lat-1):log(atc)) exp(log(atc)+(log(atc)-log(0.01))/(Lat-1):(log(atc)-log(0.01))/(Lat-1):log(atLimit)) inf]
M_R=[0 exp(log(0.01):(log(atc)-log(0.01))/(Lat-1):log(atc)) inf];

figure(1)
DDCheck(at_R);
xlabel('discretization of a','fontsize',15);
set(gca,'fontsize',15)

LalphaY=length(alphaY_R)-1
LY=length(Y_R)-1
Lat=length(at_R)-1
LM=length(M_R)-1
% Generate the Junction tree CPT and engine
[SBNengine alphaY_CPD Y1_CPD Ytc_CPD a1_CPD at_CPD M_CPD]=Crack_SBNengine_M(T,alphaY_R,Y_R,at_R,M_R,alphaY_mean,alphaY_std,Y_mean,Y_std,a0_mean,M_sigma,DeltaN,DeltaS,m,C,atc);
ev_R=FindInt(at_R,ev);
[V Iatc]=min(abs(at_R-atc));
[SBNengine_ev, loglik_ev] = enter_evidence(SBNengine, ev_R);
alphaYmarg_Jtree=zeros(LalphaY,T);
Ymarg_Jtree=zeros(LY,T);
atmarg_Jtree=zeros(Lat,T);
Mmarg_Jtree=zeros(LM,T);
for t=1:T
    marg_ev=marginal_nodes(SBNengine_ev,t*4-3);
    alphaYmarg_Jtree(:,t)=marg_ev.T;
    marg_ev=marginal_nodes(SBNengine_ev,t*4-2);
    Ymarg_Jtree(:,t)=marg_ev.T;
    marg_ev=marginal_nodes(SBNengine_ev,t*4-1);
    atmarg_Jtree(:,t)=marg_ev.T;
    marg_ev=marginal_nodes(SBNengine_ev,t*4);
    if length(marg_ev.T)==1
        marg_ev.T=zeros(LM,1);
    end
    Mmarg_Jtree(:,t)=marg_ev.T;
    FailProb_Jtree(t)=sum(atmarg_Jtree(Iatc:end,t));
end

atmarg_JtreeFine=atmarg_Jtree;
Beta_Jtree=-norminv(FailProb_Jtree);

Beta_JtreeFine=Beta_Jtree
at_RFine=at_R;
save('StochasticCrack_Fine','atmarg_JtreeFine','Beta_JtreeFine','at_RFine')