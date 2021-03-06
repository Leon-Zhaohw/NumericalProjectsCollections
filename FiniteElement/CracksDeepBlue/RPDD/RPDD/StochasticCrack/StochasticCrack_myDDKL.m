function [alphaYmarg_myDD, Ymarg_myDD, atmarg_myDD, Mmarg_myDD, alphaY_R, Y_R, at_R, M_R]=...
    StochasticCrack_myDDKL(T,DeltaN, DeltaS, m, C, atc, Inspect, alphaY_mean, alphaY_std, Y_mean, Y_std, a0_mean, M_sigma, alphaY_R, Y_R, at_R, M_R,LalphaY_Targ, LY_Targ, Lat_Targ)


% Stochastic Cratck Growth Model using Straub's DBN model with proposed Dynamic
% Discretization Scheme and KL error judgement 

disp('************** DD Algorithm with KL error ******************');
TotalIter=20;

for Iter=1:TotalIter
    LalphaY=length(alphaY_R)-1;
    LY=length(Y_R)-1;
    Lat=length(at_R)-1;
    LM=length(M_R)-1;
    % Generate the Junction tree CPT and engine
    [SBNengine alphaY_CPD Y1_CPD Ytc_CPD a1_CPD at_CPD M_CPD]=Crack_SBNengine_M(T,alphaY_R,Y_R,at_R,M_R,alphaY_mean,alphaY_std,Y_mean,Y_std,a0_mean,M_sigma,DeltaN,DeltaS,m,C,atc);
    
    %Choose slices to cal marginal prob. slices should be different with
    numnodes=4;
    ev=cell(1,T*numnodes);
    %Insert the evidence   
    InpSize=size(Inspect);
    for i=1:InpSize(1)
       ev{Inspect(i,1)*numnodes}=Inspect(i,2);
    end    
    
    ev_R=FindInt(M_R,ev);
    [V Iatc]=min(abs(at_R-atc));
    [SBNengine_ev, loglik_ev] = enter_evidence(SBNengine, ev_R);
    alphaYmarg_myDD=zeros(LalphaY,T);
    Ymarg_myDD=zeros(LY,T);
    atmarg_myDD=zeros(Lat,T);
    Mmarg_myDD=zeros(LM,T);
    for t=1:T
        marg_ev=marginal_nodes(SBNengine_ev,t*4-3);
        alphaYmarg_myDD(:,t)=marg_ev.T;
        marg_ev=marginal_nodes(SBNengine_ev,t*4-2);
        Ymarg_myDD(:,t)=marg_ev.T;
        marg_ev=marginal_nodes(SBNengine_ev,t*4-1);
        atmarg_myDD(:,t)=marg_ev.T;
        marg_ev=marginal_nodes(SBNengine_ev,t*4);
        if length(marg_ev.T)==1
            marg_ev.T=zeros(LM,1);
        end
        Mmarg_myDD(:,t)=marg_ev.T;
    end
    
    if (LalphaY>=LalphaY_Targ && LY>=LY_Targ && Lat>=Lat_Targ) || (Iter>2 && isempty(alphaYAdd) && isempty(YAdd) && isempty(atAdd))
        break
    end
    
    % Insert Virtual Inspect
    atmarg2_myDD=cell(1);
    Ymarg2_myDD=cell(1);
    alphaYmarg2_myDD=cell(1);
    Count1=1;
    for t=2:max(Inspect(:,1))
        ev2=ev_R;
        for VirtualEv=Lat:-1:1
             if atmarg_myDD(VirtualEv,t)>0
                 break
             end
        end
        ev2{1,t*4-1}=VirtualEv;
        [engine_ev,loglik_ev]=enter_evidence(SBNengine,ev2);
        Count2=1;
        Count3=1;
        for tt=1:t-1   
            marg_ev=marginal_nodes(engine_ev,tt*4-1);
            atmarg2_myDD{Count1}(:,Count2)=marg_ev.T; 
            marg_ev=marginal_nodes(engine_ev,tt*4-2);
            Ymarg2_myDD{Count1}(:,Count2)=marg_ev.T;            
            Count2=Count2+1;
        end
        for tt=1:max(Inspect(:,1))
            marg_ev=marginal_nodes(engine_ev,tt*4-3);
            alphaYmarg2_myDD{Count1}(:,Count3)=marg_ev.T;
            Count3=Count3+1;
        end
        Count1=Count1+1;
    end

    %%%%%%%%%%%%%%%%%%%%%% partition step %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    alphaYMerg=[];
    YMerg=[];
    atMerg=[];
    MMerg=[];
    alphaYAdd=[];
    YAdd=[];
    atAdd=[];
    MMerg=[];

    %disp('%%%%%%%%%%%%%%%%%%%%%%%% Add Step %%%%%%%%%%%%%%%%%%%%%%%%%')
    alphaYAdd1=AddStep_KL(alphaYmarg_myDD, alphaY_R);
    alphaYAdd2=AddStep_KL(alphaYmarg2_myDD, alphaY_R);
    YAdd1=AddStep_KL(Ymarg_myDD,Y_R);
    YAdd2=AddStep_KL(Ymarg2_myDD,Y_R);
    atAdd1=AddStep_KL(atmarg_myDD,at_R);
    disp('******************************************************************************')
    atAdd2=AddStep_KL(atmarg2_myDD,at_R); 
    disp('*******************************************************************************')
    alphaYAdd=[alphaYAdd1 alphaYAdd2];
    alphaYAdd=sort(unique(alphaYAdd),'ascend');
    if alphaYAdd(end)==LalphaY
        alphaYAdd(end)=[];
    end       
    if length(alphaYAdd)+length(alphaY_R)-1>=LalphaY_Targ
        alphaYAdd=randsample(alphaYAdd,LalphaY_Targ-length(alphaY_R)+1);
    end       

    YAdd=[YAdd1 YAdd2];
    YAdd=sort(unique(YAdd),'ascend');
    if YAdd(end)==LY
        YAdd(end)=[];
    end          
    if length(YAdd)+length(Y_R)-1>=LY_Targ
        YAdd=randsample(YAdd,LY_Targ-length(Y_R)+1);
    end

    atAdd=[atAdd1 atAdd2];
    atAdd=sort(unique(atAdd),'ascend');
    
    if atAdd(end)==Lat
        atAdd(end)=[];
    end
    if length(atAdd)+length(at_R)-1>=Lat_Targ
        atAdd=randsample(atAdd,Lat_Targ-length(at_R)+1);
    end
    if isempty(atAdd)
        MAdd=[];
    else
        MAdd=atAdd;
        for i=length(MAdd):-1:1
            if MAdd(i)>Iatc
                MAdd(i)=[];
            end            
        end
        if isempty(MAdd) || MAdd(end)==LM
            MAdd(end)=[];
        end 
    end
        
    % Update the discretization
    if isempty(alphaYAdd)
        alphaYAdd=[];
    end
    if isempty(YAdd)
        YAdd=[];
    end
    if isempty(atAdd)
        atAdd=[];
    end
    if isempty(MAdd)
        MAdd=[];
    end
    
    alphaY_R=DDStep(alphaYAdd,alphaYMerg,alphaY_R);
    Y_R=DDStep(YAdd, YMerg,Y_R);
    at_R=DDStep(atAdd, atMerg,at_R);
    M_R=DDStep(MAdd,MMerg,M_R);
    
    
end
 













