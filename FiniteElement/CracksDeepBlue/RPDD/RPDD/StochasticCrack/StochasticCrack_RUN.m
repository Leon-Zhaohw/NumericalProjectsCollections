clear all;
clc;
format long e
T=50;

%Choose if want to add the Measurement nodes(1) or not(2)
MeasureNode=1;

%Inspection DATA[slice Cratck Length]
Inspect=[10 1.4
        20 2
        30 5];

DeltaN=1E4; 
DeltaS=60;
m=3.59;
C=1.42E-8/1000^(3.59/2);
atc=60;
atLimit=99;
%alphaY node
alphaY_mean=1E4;
alphaY_std=1E4;
%Y node
Y_mean=1.2;
Y_std=0.8;
%a node
a0_mean=1;
%M node
M_sigma=5;



%Choose slices to cal marginal prob. slices should be different with
numnodes=4;
ev=cell(1,T*numnodes);
%Insert the evidence   
InpSize=size(Inspect);
for i=1:InpSize(1)
   ev{Inspect(i,1)*numnodes}=Inspect(i,2);
end

% run the Junction tree model with very dense intervals
%[atmarg_JtreeFine Beta_JtreeFine,at_RFine,Ytc_CPD,at_CPD,M_CPD]=StochasticCrack_Fine(alphaY_mean,alphaY_std,Y_mean,Y_std,a0_mean,M_sigma,DeltaN,DeltaS,m,C,atc,atLimit,ev,T);
%break
load StochasticCrack_Fine


LalphaY=4;
LY=4;
Lat=4;
LM=Lat;
alphaY_R_Init=[0 10.^(1:(6-1)/(LalphaY-1):6)];
Y_R_Init=[0.0:15/(LY-1):15 inf];
at_R_Init=[0 exp(log(0.01):(log(atc)-log(0.01))/(Lat-1):log(atc)) atLimit inf];
%at_R=[0:atLimit/(Lat-1):atLimit inf];
[V Iatc]=min(abs(at_R_Init-atc));
at_R_Init(Iatc(1))=atc; 
M_R_Init=[at_R_Init(1:end-2) inf];

alphaY_R_myDD=alphaY_R_Init;
Y_R_myDD=Y_R_Init;
at_R_myDD=at_R_Init;
M_R_myDD=M_R_Init;
alphaY_R_myDDKL=alphaY_R_Init;
Y_R_myDDKL=Y_R_Init;
at_R_myDDKL=at_R_Init;
M_R_myDDKL=M_R_Init;

TotalIter=9;
 
for Iter=1:TotalIter
    Iter
    LalphaY(Iter)=length(alphaY_R_myDD)-1;
    LY(Iter)=length(Y_R_myDD)-1;
    Lat(Iter)=length(at_R_myDD)-1;
    LM(Iter)=length(M_R_myDD)-1;    

    [alphaYmarg_myDD, Ymarg_myDD, atmarg_myDD, Mmarg_myDD, alphaY_R_myDD2, Y_R_myDD2, at_R_myDD2, M_R_myDD2]=...
        StochasticCrack_myDD(T,DeltaN, DeltaS, m, C, atc,Inspect, alphaY_mean, alphaY_std, Y_mean, Y_std, a0_mean, M_sigma, alphaY_R_myDD, Y_R_myDD, at_R_myDD, M_R_myDD);
     
    [alphaYmarg_Unif,Ymarg_Unif, atmarg_Unif, Mmarg_Unif,alphaY_R_Unif, Y_R_Unif, at_R_Unif,M_R_Unif,alphaYmarg_Log,Ymarg_Log,atmarg_Log, Mmarg_Log, alphaY_R_Log, Y_R_Log, at_R_Log, M_R_Log]=...
        StochasticCrack_Unif_Log(T,DeltaN, DeltaS, m, C, atc,atLimit, Inspect, alphaY_mean, alphaY_std, Y_mean, Y_std, a0_mean, M_sigma, LalphaY(Iter), LY(Iter), Lat(Iter));
     
    [alphaYmarg_myDDKL, Ymarg_myDDKL, atmarg_myDDKL, Mmarg_myDDKL, alphaY_R_myDDKL2, Y_R_myDDKL2, at_R_myDDKL2, M_R_myDDKL2]=...
    StochasticCrack_myDDKL(T,DeltaN, DeltaS, m, C, atc, Inspect, alphaY_mean, alphaY_std, Y_mean, Y_std, a0_mean, M_sigma, alphaY_R_myDDKL, Y_R_myDDKL, at_R_myDDKL, M_R_myDDKL,LalphaY(Iter), LY(Iter), Lat(Iter));
     
    alphaYmarg_Cell{Iter}={alphaYmarg_myDD,alphaYmarg_myDDKL,alphaYmarg_Unif,alphaYmarg_Log};
    Ymarg_Cell{Iter}={Ymarg_myDD,Ymarg_myDDKL,Ymarg_Unif,Ymarg_Log};
    atmarg_Cell{Iter}={atmarg_myDD,atmarg_myDDKL, atmarg_Unif,atmarg_Log};
    Mmarg_Cell{Iter}={Mmarg_myDD,Mmarg_myDDKL, Mmarg_Unif, Mmarg_Log};
    alphaY_R_Cell{Iter}={alphaY_R_myDD,alphaY_R_myDDKL2,alphaY_R_Unif,alphaY_R_Log};
    Y_R_Cell{Iter}={Y_R_myDD,Y_R_myDDKL2,Y_R_Unif,Y_R_Log};
    at_R_Cell{Iter}={at_R_myDD,at_R_myDDKL2, at_R_Unif,at_R_Log};
    M_R_Cell{Iter}={M_R_myDD,M_R_myDDKL2, M_R_Unif, M_R_Log};      
    alphaY_R_myDD=alphaY_R_myDD2;
    Y_R_myDD=Y_R_myDD2;
    at_R_myDD=at_R_myDD2;
    M_R_myDD=M_R_myDD2;
    alphaY_R_myDDKL=alphaY_R_myDDKL2;
    Y_R_myDDKL=Y_R_myDDKL2;
    at_R_myDDKL=at_R_myDDKL2;
    M_R_myDDKL=M_R_myDDKL2;    

end
save('Output','alphaYmarg_Cell','Ymarg_Cell','atmarg_Cell','Mmarg_Cell','alphaY_R_Cell','Y_R_Cell','at_R_Cell','M_R_Cell','Lat');
load('Output','alphaYmarg_Cell','Ymarg_Cell','atmarg_Cell','Mmarg_Cell','alphaY_R_Cell','Y_R_Cell','at_R_Cell','M_R_Cell','Lat');



for Iter=1:TotalIter
    %% Beta Error
    for t=1:T
        Beta_myDD(t)=ReliabIndex(atmarg_Cell{Iter}{1}(:,t),at_R_Cell{Iter}{1},atc);
        Beta_myDDKL(t)=ReliabIndex(atmarg_Cell{Iter}{2}(:,t),at_R_Cell{Iter}{2},atc);
        Beta_Unif(t)=ReliabIndex(atmarg_Cell{Iter}{3}(:,t),at_R_Cell{Iter}{3},atc);
        Beta_Log(t)=ReliabIndex(atmarg_Cell{Iter}{4}(:,t),at_R_Cell{Iter}{4},atc);
        
    end
    BetaError_myDD(Iter)=sum(abs(Beta_JtreeFine-Beta_myDD))/T;
    BetaError_myDDKL(Iter)=sum(abs(Beta_JtreeFine-Beta_myDDKL))/T;
    BetaError_Unif(Iter)=sum(abs(Beta_JtreeFine-Beta_Unif))/T;
    BetaError_Log(Iter)=sum(abs(Beta_JtreeFine-Beta_Log))/T;
    
    %% KL Error
    for t=1:T
        klerror_myDD(t)=sum(KLError(at_R_Cell{Iter}{1},atmarg_Cell{Iter}{1}(:,t)));
        klerror_myDDKL(t)=sum(KLError(at_R_Cell{Iter}{2},atmarg_Cell{Iter}{2}(:,t)));
        klerror_Unif(t)=sum(KLError(at_R_Cell{Iter}{3},atmarg_Cell{Iter}{3}(:,t)));  
        klerror_Log(t)=sum(KLError(at_R_Cell{Iter}{4},atmarg_Cell{Iter}{4}(:,t))); 
       
    end
    KLError_myDD(Iter)=mean(klerror_myDD);
    KLError_myDDKL(Iter)=mean(klerror_myDDKL);
    KLError_Unif(Iter)=mean(klerror_Unif);
    KLError_Log(Iter)=mean(klerror_Log);  
end
% Initial Crack size inference
[a1_R_Fine,a1cdf_Fine]=InitCrack(at_RFine,atmarg_JtreeFine(:,1));
[a1_R_myDD,a1cdf_myDD]=InitCrack(at_R_Cell{TotalIter}{1},atmarg_Cell{TotalIter}{1}(:,1));
[a1_R_myDDKL,a1cdf_myDDKL]=InitCrack(at_R_Cell{TotalIter}{2},atmarg_Cell{TotalIter}{2}(:,1));
[a1_R_Unif,a1cdf_Unif]=InitCrack(at_R_Cell{TotalIter}{3},atmarg_Cell{TotalIter}{3}(:,1));
[a1_R_Log,a1cdf_Log]=InitCrack(at_R_Cell{TotalIter}{4},atmarg_Cell{TotalIter}{4}(:,1));

% Plot final discretization
h1=figure('visible','off');
subplot(4,1,1);
DDCheck(alphaY_R_Cell{end}{1});
xlabel('Discretization of \alpha_X_t','fontsize',15);
set(gca,'fontsize',15,'ytick',[]);
subplot(4,1,2);
DDCheck(Y_R_Cell{end}{1});
xlabel('Discretization of X_t','fontsize',15);
set(gca,'fontsize',15,'ytick',[]);
subplot(4,1,3);
DDCheck(at_R_Cell{end}{1});
xlabel('Discretization of A_t','fontsize',15);
set(gca,'fontsize',15,'ytick',[]);
subplot(4,1,4);
DDCheck(M_R_Cell{end}{1});
xlabel('Discretization of M_t','fontsize',15);
set(gca,'fontsize',15,'ytick',[]);
saveas(h1,'figures/StochasticCrackDiscrete','fig');
print(h1,'-depsc','figures/StochasticCrackDiscrete.eps');

%% Plot the Beta Error
h2=figure('visible','off'); hold on;
plot(Lat(1:TotalIter),log10(BetaError_myDD),'--sr','linewidth',2);
plot(Lat(1:TotalIter),log10(BetaError_myDDKL),'--og','linewidth',2);
plot(Lat(1:TotalIter),log10(BetaError_Unif),'--<b','linewidth',2);
plot(Lat(1:TotalIter),log10(BetaError_Log),'--*m','linewidth',2);
legend('RPDD','RPDD with KL error','Uniform','Logarithm');
xlabel('Number of intervals for A_t','fontsize',15)
ylabel('log_1_0(\epsilon_\beta)','fontsize',15);
set(gca,'fontsize',15)
grid on 
saveas(h2,'figures/StochasticCrackBetaError','fig');
print(h2,'-depsc','figures/StochasticCrackBetaError.eps');


%% Plot the KL distance Error
h3=figure('visible','off'); hold on;
plot(Lat(1:TotalIter),log10(KLError_myDD),'--sr','linewidth',2);
plot(Lat(1:TotalIter),log10(KLError_myDDKL),'--og','linewidth',2);
plot(Lat(1:TotalIter),log10(KLError_Unif),'--<b','linewidth',2);
plot(Lat(1:TotalIter),log10(KLError_Log),'--*m','linewidth',2);

legend('RPDD','RPDD with KL error','Uniform','Logarithm');
xlabel('Number of intervals for A_t','fontsize',15)
ylabel('log_1_0(\epsilon_K_L)','fontsize',15);
set(gca,'fontsize',15)
grid on 
saveas(h3,'figures/StochasticCrackKLError','fig');
print(h3,'-depsc','figures/StochasticCrackKLError.eps');


%% Plot the failure probability
h4=figure('visible','off'); hold on;
plot(1:T,Beta_myDD,'--r','linewidth',2);
plot(1:T,Beta_myDDKL,'--g','linewidth',2);
plot(1:T,Beta_Unif,'--b','linewidth',2);
plot(1:T,Beta_Log,'--m','linewidth',2);
plot(1:T,Beta_JtreeFine,'-k','linewidth',2);
legend('RPDD ','RPDD with KL error','Uniform','Logarithm','Logarithm, dense');
xlabel('Time step','fontsize',15)
ylabel('\beta','fontsize',15)
set(gca,'fontsize',15);
grid on
saveas(h4,'figures/StochasticCrackBeta','fig');
print(h4,'-depsc','figures/StochasticCrackBeta.eps');


%plot the cdf distribution initial crack size
h5=figure('visible','off'); hold on;
plot(a1_R_myDD,a1cdf_myDD,'--r','linewidth',2);
plot(a1_R_myDDKL,a1cdf_myDDKL,'--g','linewidth',2);
plot(a1_R_Unif,a1cdf_Unif,'--b','linewidth',2);
plot(a1_R_Log,a1cdf_Log,'--m','linewidth',2);
plot(a1_R_Fine,a1cdf_Fine,'-k','linewidth',2);
hold on
xlabel('initial crack size','fontsize',15);
ylabel('CDF','fontsize',15);
legend('RPDD ','RPDD with KL error','Uniform','Logarithm','Logarithm, dense');
set(gca,'fontsize',15)
axis([0 3 0 1])
grid on
saveas(h5,'figures/StochasticCrackInitialCrack','fig');
print(h5,'-depsc','figures/StochasticCrackInitialCrack.eps');
























