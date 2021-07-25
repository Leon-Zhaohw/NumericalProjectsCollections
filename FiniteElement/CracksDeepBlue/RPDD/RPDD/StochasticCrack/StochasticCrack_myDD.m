function [alphaYmarg_myDD, Ymarg_myDD, atmarg_myDD, Mmarg_myDD, alphaY_R, Y_R, at_R, M_R]=...
    StochasticCrack_myDD(T,DeltaN, DeltaS, m, C, atc, Inspect, alphaY_mean, alphaY_std, Y_mean, Y_std, a0_mean, M_sigma, alphaY_R, Y_R, at_R, M_R)


% Stochastic Cratck Growth Model using Straub's DBN model with Dynamic
% Discretization Scheme



MAdd=[];
MMerg=[];
atAdd=[];
atMerg=[];
YAdd=[];
YMerg=[];
alphaYAdd=[];
alphaYMerg=[];

%Choose slices to cal marginal prob. slices should be different with
numnodes=4;
ev=cell(1,T*numnodes);
%Insert the evidence   
InpSize=size(Inspect);
for i=1:InpSize(1)
   ev{Inspect(i,1)*numnodes}=Inspect(i,2);
end

disp('********************* DD Algorithm **********************');

LalphaY=length(alphaY_R)-1;
LY=length(Y_R)-1;
Lat=length(at_R)-1;
LM=length(M_R)-1;



% Generate the Junction tree CPT and engine
[SBNengine alphaY_CPD Y1_CPD Ytc_CPD a1_CPD at_CPD M_CPD]=Crack_SBNengine_M(T,alphaY_R,Y_R,at_R,M_R,alphaY_mean,alphaY_std,Y_mean,Y_std,a0_mean,M_sigma,DeltaN,DeltaS,m,C,atc);

ev_R=FindInt(at_R,ev);
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
[alphaYAdd_Smooth alphaYE_Smooth]=AddStep_Smooth(alphaYmarg2_myDD);
[YAdd_Smooth YE_Smooth]=AddStep_Smooth(Ymarg2_myDD);
[atAdd_Smooth atE_Smooth]=AddStep_Smooth(atmarg2_myDD);
[YAdd_Predict YE_Predict atAdd_Predict atE_Predict]=AddStep_Predict(atmarg_myDD,Ymarg_myDD,at_CPD,T,Inspect);
atAdd_Smooth=sort(unique(atAdd_Smooth),'ascend');
atAdd_Predict=sort(unique(atAdd_Predict),'ascend');
atAdd=[atAdd_Predict atAdd_Smooth];      
YAdd_Smooth=sort(unique(YAdd_Smooth),'ascend');
YAdd_Predict=sort(unique(YAdd_Predict),'ascend');       
YAdd=[YAdd_Predict YAdd_Smooth]; 
alphaYAdd=[alphaYAdd_Smooth];
alphaYAdd=sort(unique(alphaYAdd),'ascend');
if alphaYAdd(end)==LalphaY
    alphaYAdd(end)=[];
end  
YAdd=sort(unique(YAdd),'ascend');
if YAdd(end)==LY
    YAdd(end)=[];
end  
atAdd=sort(unique(atAdd),'ascend');
if atAdd(end)==Lat
    atAdd(end)=[];
end     
MAdd=atAdd;
for i=length(MAdd):-1:1
    if MAdd(i)>Iatc
        MAdd(i)=[];
    end            
end
if MAdd(end)==LM
    MAdd(end)=[];
end   

% Update the discretization
alphaY_R=DDStep(alphaYAdd,alphaYMerg,alphaY_R);
Y_R=DDStep(YAdd, YMerg,Y_R);
at_R=DDStep(atAdd, atMerg,at_R);
M_R=DDStep(MAdd,MMerg,M_R);


% 
% % Plot final distribution
% figure(6)
% subplot(4,1,1);
% DDCheck(alphaY_R);
% xlabel('discretization of \alpha_X_t','fontsize',15);
% set(gca,'fontsize',15,'ytick',[]);
% subplot(4,1,2);
% DDCheck(Y_R);
% xlabel('discretization of X_t','fontsize',15);
% set(gca,'fontsize',15,'ytick',[]);
% subplot(4,1,3);
% DDCheck(at_R);
% xlabel('discretization of a_t','fontsize',15);
% set(gca,'fontsize',15,'ytick',[]);
% subplot(4,1,4);
% DDCheck(M_R);
% xlabel('discretization of M_t','fontsize',15);
% set(gca,'fontsize',15,'ytick',[]);
% 
% figure(7)
% plot(1:T, Beta_myDDFine,'-b','linewidth',3)
% hold on
% plot(1:T, Beta_myDDRPDD(:,Iter),'--r','linewidth',3)
% xlabel('time step','fontsize',15)
% ylabel('\beta','fontsize',15)
% l1=legend(['logarithm, with ' num2str(length(at_RFine)-1) ' intervals of a_t'],['RPDD, with ' num2str(Lat(end)-1) ' intervals for a_t'])
% set(l1,'fontsize',12)
% set(gca,'fontsize',15)
% grid on
% 
% figure(8)
% %plot the cdf distribution initial crack size
% a1cdf_Fine=zeros(1,length(atmarg_myDDFine(:,1)));
% a1cdf_Fine(1)=atmarg_myDDFine(1,1);
% for i=2:length(a1cdf_Fine)
%     a1cdf_Fine(i)=a1cdf_Fine(i-1)+atmarg_myDDFine(i,1);
% end
% a1cdf=zeros(1,length(atmarg_myDD(:,1)));
% a1cdf(1)=atmarg_myDD(1,1);
% for i=2:length(a1cdf)
%     a1cdf(i)=a1cdf(i-1)+atmarg_myDD(i,1);
% end
% plot(at_RFine(1:end-1),a1cdf_Fine,'b','linewidth',3);
% hold on
% plot(at_R(1:end-1),a1cdf,'--r','linewidth',3);
% xlabel('initial crack size','fontsize',15);
% ylabel('CDF','fontsize',15);
% l2=legend(['logarithm, with ' num2str(length(at_RFine)-1) ' intervals for a_t'],['RPDD, with ' num2str(Lat(end)-1) ' intervals for a_t'])
% set(l2,'fontsize',12)
% set(gca,'fontsize',15)
% axis([0 3 0 1])
% grid on









