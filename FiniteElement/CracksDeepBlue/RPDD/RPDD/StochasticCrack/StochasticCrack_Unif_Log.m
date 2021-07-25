% This is a comparision plot for constant cratck with uniform and log discretization
% This code should be run after the "StochasticCratck_myDD.m"


function [alphaYmarg_Unif,Ymarg_Unif, atmarg_Unif, Mmarg_Unif,alphaY_R_Unif, Y_R_Unif, at_R_Unif,M_R_Unif,alphaYmarg_Log,Ymarg_Log,atmarg_Log, Mmarg_Log, alphaY_R_Log, Y_R_Log, at_R_Log, M_R_Log]=...
    StochasticCrack_Unif_Log(T,DeltaN, DeltaS, m, C, atc,atLimit, Inspect, alphaY_mean, alphaY_std, Y_mean, Y_std, a0_mean, M_sigma, LalphaY, LY, Lat)


%Choose slices to cal marginal prob. slices should be different with
numnodes=4;
ev=cell(1,T*numnodes);
%Insert the evidence   
InpSize=size(Inspect);
for i=1:InpSize(1)
   ev{Inspect(i,1)*numnodes}=Inspect(i,2);
end

disp('**** Uniform Discretization ****');
alphaY_R_Unif=[0 10^1:(10^6-10^1)/(LalphaY-1):10^6];
Y_R_Unif=[0.0:15/(LY-1):15 inf];
at_R_Unif=[0:atLimit/(Lat-1):atLimit inf];
[V Iatc]=min(abs(at_R_Unif-atc));
at_R_Unif(Iatc(1))=atc; 
M_R_Unif=[at_R_Unif(1:Iatc) inf];
LM=length(M_R_Unif)-1;
[SBNengine_Unif, alphaY_CPD_Unif, Y1_CPD_Unif, Ytc_CPD_Unif, a1_CPD_Unif, at_CPD_Unif, M_CPD_Unif]=Crack_SBNengine_M(T,alphaY_R_Unif,Y_R_Unif,at_R_Unif,M_R_Unif,alphaY_mean,alphaY_std,Y_mean,Y_std,a0_mean,M_sigma,DeltaN,DeltaS,m,C,atc);
ev_R=FindInt(M_R_Unif,ev);
[SBNengine_ev, loglik_ev] = enter_evidence(SBNengine_Unif, ev_R);
alphaYmarg_Unif=zeros(LalphaY,T);
Ymarg_Unif=zeros(LY,T);
atmarg_Unif=zeros(Lat,T);
Mmarg_Unif=zeros(LM,T);
for t=1:T
    marg_ev=marginal_nodes(SBNengine_ev,t*4-3);
    alphaYmarg_Unif(:,t)=marg_ev.T;
    marg_ev=marginal_nodes(SBNengine_ev,t*4-2);
    Ymarg_Unif(:,t)=marg_ev.T;
    marg_ev=marginal_nodes(SBNengine_ev,t*4-1);
    atmarg_Unif(:,t)=marg_ev.T;
    marg_ev=marginal_nodes(SBNengine_ev,t*4);
    if length(marg_ev.T)==1
        marg_ev.T=zeros(LM,1);
    end
    Mmarg_Unif(:,t)=marg_ev.T;
end


disp('**** Logrithmic Discretization ****');
alphaY_R_Log=[0 10.^(1:5/(LalphaY-1):6)];
Y_R_Log=[0.0:15/(LY-1):15 inf];
at_R_Log=[0 exp(log(0.01):(log(atLimit)-log(0.01))/(Lat-2):log(atLimit)) inf];
[V Iatc]=min(abs(at_R_Log(1:end-2)-atc));
at_R_Log(Iatc(1))=atc;   
M_R_Log=[at_R_Log(1:Iatc) inf];
LM=length(M_R_Log)-1;
[SBNengine_Log, alphaY_CPD_Log, Y1_CPD_Log, Ytc_CPD_Log, a1_CPD_Log, at_CPD_Log, M_CPD_Log,]=Crack_SBNengine_M(T,alphaY_R_Log,Y_R_Log,at_R_Log,M_R_Log,alphaY_mean,alphaY_std,Y_mean,Y_std,a0_mean,M_sigma,DeltaN,DeltaS,m,C,atc);
ev_R=FindInt(M_R_Log,ev);
[SBNengine_ev, loglik_ev] = enter_evidence(SBNengine_Log, ev_R);
alphaYmarg_Log=zeros(LalphaY,T);
Ymarg_Log=zeros(LY,T);
atmarg_Log=zeros(Lat,T);
Mmarg_Log=zeros(LM,T);
for t=1:T
    marg_ev=marginal_nodes(SBNengine_ev,t*4-3);
    alphaYmarg_Log(:,t)=marg_ev.T;
    marg_ev=marginal_nodes(SBNengine_ev,t*4-2);
    Ymarg_Log(:,t)=marg_ev.T;
    marg_ev=marginal_nodes(SBNengine_ev,t*4-1);
    atmarg_Log(:,t)=marg_ev.T;
    marg_ev=marginal_nodes(SBNengine_ev,t*4);
    if length(marg_ev.T)==1
        marg_ev.T=zeros(LM,1);
    end
    Mmarg_Log(:,t)=marg_ev.T;
end
