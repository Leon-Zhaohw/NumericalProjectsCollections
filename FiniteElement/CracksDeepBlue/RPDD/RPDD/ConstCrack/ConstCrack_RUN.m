clc
clear all
format long
% T is the time slice
T=20;
% mean and std value of initial crak size
mu_a0=7;
sigma_a0=3;
% mean and std of da
mu_da=1;
sigma_da=0.3;
% measurment error
sigma_error=2;
atc=30;

% Inspection data
% [t1 Cracksize1]
 Evidence=[5 10;
           10 14
           15 25];

% discretization range for at and da

daLimit=5;
atLimit=50;
Lda=4;
Lat=4;
LM=4;

TotalIter=9;

% Initial discretization
% the discretization limit for at, different with atc which is the
% critical crack size;
da_R_Init=[0:daLimit/(Lda-1):daLimit inf];
at_R_Init=[0:atc/(Lat-1):atc atc+atc/(Lat-1):atc/(Lat-1):atLimit inf];
M_R_Init=[0:atc/(Lat-1):atc inf];

% Bayesian Linear Regression Method
load FailProb_BLR
 %[FailProb_BLR mu_at_BLR sigma2_at_BLR mu_da_BLR sigma2_da_BLR]=ConstCrack_BLR(Evidence,T,mu_a0,sigma_a0,mu_da,sigma_da,sigma_error,atc)
Beta_BLR=-norminv(FailProb_BLR);
% Kalman Smoother 
load FailProb_KS
%[FailProb_KS mu_at_KS sigma2_at_KS]=ConstCrack_KF(Evidence,T,mu_a0,sigma_a0,mu_da,sigma_da,sigma_error,atc)

da_R_myDD=da_R_Init;
at_R_myDD=at_R_Init;
M_R_myDD=M_R_Init;
da_R_myDDKL=da_R_Init;
at_R_myDDKL=at_R_Init;
M_R_myDDKL=M_R_Init;
for t=1:T
    da_R_DDNeil{t}=da_R_Init;
    at_R_DDNeil{t}=at_R_Init;
    M_R_DDNeil{t}=M_R_Init;
end
a0_R_DDNeil=at_R_Init;

for Iter=1:TotalIter
    Iter    
    Lda(Iter)=length(da_R_myDD)-1;
    Lat(Iter)=length(at_R_myDD)-1;
    LM(Iter)=length(M_R_myDD)-1;
    disp('**** RPDD ****');
    time0=cputime;
    [damarg_myDD, atmarg_myDD, Mmarg_myDD, da_R_myDD2, at_R_myDD2, M_R_myDD2]=ConstCrack_myDD (mu_a0, sigma_a0, mu_da, sigma_da, sigma_error,Evidence, atc,T, da_R_myDD,at_R_myDD,M_R_myDD);
    time_RPDD(Iter)=cputime-time0
    disp('**** Uniform ****');
    da_R_Unif=[0:daLimit/(Lda(Iter)-1):daLimit inf];
    at_R_Unif=[0:atLimit/(Lat(Iter)-1):atLimit inf];
    [V Iatc]=min(abs(at_R_Unif-atc));
    at_R_Unif(Iatc(1))=atc;     
    M_R_Unif=[at_R_Unif(1:Iatc(1)) inf];
    time0=cputime;
    [damarg_Unif, atmarg_Unif, Mmarg_Unif]=ConstCrack_Unif(mu_a0, sigma_a0, mu_da, sigma_da, sigma_error,Evidence, atc, T,da_R_Unif, at_R_Unif,M_R_Unif);
    time_Unif(Iter)=cputime-time0
    disp('**** Neil Discretization ****');
    [a0marg_DDNeil, damarg_DDNeil, atmarg_DDNeil, Mmarg_DDNeil, da_R_DDNeil, at_R_DDNeil, M_R_DDNeil]...
        =CostCrack_DDNeil(mu_a0, sigma_a0, mu_da, sigma_da, sigma_error,Evidence, atc, T, da_R_DDNeil,at_R_DDNeil,M_R_DDNeil,a0_R_DDNeil,Lda(Iter),Lat(Iter), LM(Iter));
    disp('****RPDD with KL Error ****');
    [damarg_myDDKL, atmarg_myDDKL, Mmarg_myDDKL, da_R_myDDKL2, at_R_myDDKL2, M_R_myDDKL2]...
        =ConstCrack_myDDKL(mu_a0, sigma_a0, mu_da, sigma_da, sigma_error,Evidence, atc,T, da_R_myDDKL,at_R_myDDKL,M_R_myDDKL,Lda(Iter), Lat(Iter));
    
    damarg_Cell{Iter}={damarg_myDD,damarg_Unif,damarg_DDNeil,damarg_myDDKL};
    atmarg_Cell{Iter}={atmarg_myDD,atmarg_Unif,atmarg_DDNeil,atmarg_myDDKL};
    Mmarg_Cell{Iter}={Mmarg_myDD,Mmarg_Unif, Mmarg_DDNeil,Mmarg_myDDKL};
    da_R_Cell{Iter}={da_R_myDD,da_R_Unif,da_R_DDNeil, da_R_myDDKL2};
    at_R_Cell{Iter}={at_R_myDD,at_R_Unif,at_R_DDNeil,at_R_myDDKL2};
    M_R_Cell{Iter}={M_R_myDD,M_R_Unif,M_R_DDNeil, M_R_myDDKL2}; 
    

    da_R_myDD=da_R_myDD2;
    at_R_myDD=at_R_myDD2;
    M_R_myDD=M_R_myDD2;
    da_R_myDDKL=da_R_myDDKL2;
    at_R_myDDKL=at_R_myDDKL2;
    M_R_myDDKL=M_R_myDDKL2;      
end

%save('Output','damarg_Cell','atmarg_Cell','Mmarg_Cell','da_R_Cell','at_R_Cell','M_R_Cell','Lat');
load('Output','damarg_Cell','atmarg_Cell','Mmarg_Cell','da_R_Cell','at_R_Cell','M_R_Cell','Lat');

TotalIter=9
for Iter=1:TotalIter
    %% Beta Error
    for t=1:T
        Beta_myDD(t)=ReliabIndex(atmarg_Cell{Iter}{1}(:,t),at_R_Cell{Iter}{1},atc);
        Beta_Unif(t)=ReliabIndex(atmarg_Cell{Iter}{2}(:,t),at_R_Cell{Iter}{2},atc);
        Beta_DDNeil(t)=ReliabIndex(atmarg_Cell{Iter}{3}{t},at_R_Cell{Iter}{3}{t},atc);    
        Beta_myDDKL(t)=ReliabIndex(atmarg_Cell{Iter}{4}(:,t),at_R_Cell{Iter}{4},atc);
    end
    BetaError_myDD(Iter)=sum(abs(Beta_BLR-Beta_myDD))/T;
    BetaError_Unif(Iter)=sum(abs(Beta_BLR-Beta_Unif))/T;
    BetaError_DDNeil(Iter)=sum(abs(Beta_BLR-Beta_DDNeil))/T;
    BetaError_myDDKL(Iter)=sum(abs(Beta_BLR-Beta_myDDKL))/T;
    %% KL Error
    for t=1:T
        klerror_myDD(t)=sum(KLError(at_R_Cell{Iter}{1},atmarg_Cell{Iter}{1}(:,t)));
        klerror_Unif(t)=sum(KLError(at_R_Cell{Iter}{2},atmarg_Cell{Iter}{2}(:,t)));  
        klerror_DDNeil(t)=sum(KLError(at_R_Cell{Iter}{3}{t},atmarg_Cell{Iter}{3}{t}));  
        klerror_myDDKL(t)=sum(KLError(at_R_Cell{Iter}{4},atmarg_Cell{Iter}{4}(:,t)));
    end
    KLError_myDD(Iter)=mean(klerror_myDD);
    KLError_Unif(Iter)=mean(klerror_Unif);
    KLError_DDNeil(Iter)=mean(klerror_DDNeil);
    KLError_myDDKL(Iter)=mean(klerror_myDDKL);
    % Mean and Std 
    for t=1:T
       [mu_da_myDD(t),std_da_myDD(t)]=MeanStd(damarg_Cell{Iter}{1}(:,t),da_R_Cell{Iter}{1});
       [mu_da_Unif(t),std_da_Unif(t)]=MeanStd(damarg_Cell{Iter}{2}(:,t),da_R_Cell{Iter}{2});
       [mu_da_DDNeil(t),std_da_DDNeil(t)]=MeanStd(damarg_Cell{Iter}{3}{t},da_R_Cell{Iter}{3}{t});
       [mu_da_myDDKL(t),std_da_myDDKL(t)]=MeanStd(damarg_Cell{Iter}{4}(:,t),da_R_Cell{Iter}{4});
       [mu_at_myDD(t),std_at_myDD(t)]=MeanStd(atmarg_Cell{Iter}{1}(:,t),at_R_Cell{Iter}{1});
       [mu_at_Unif(t),std_at_Unif(t)]=MeanStd(atmarg_Cell{Iter}{2}(:,t),at_R_Cell{Iter}{2});
       [mu_at_DDNeil(t),std_at_DDNeil(t)]=MeanStd(atmarg_Cell{Iter}{3}{t},at_R_Cell{Iter}{3}{t});    
       [mu_at_myDDKL(t),std_at_myDDKL(t)]=MeanStd(atmarg_Cell{Iter}{4}(:,t),at_R_Cell{Iter}{4});
    end    
end


%% Plot the discretization
h1=figure('visible','off');
subplot(3,1,1);
DDCheck(da_R_Cell{TotalIter}{1});
xlabel('Discretization of dA_t (mm)','fontsize',15);
set(gca,'fontsize',15,'ytick',[]);
subplot(3,1,2);
DDCheck(at_R_Cell{TotalIter}{1});
xlabel('Discretization of A_t (mm)','fontsize',15);
set(gca,'fontsize',15,'ytick',[]);
subplot(3,1,3);
DDCheck(M_R_Cell{TotalIter}{1});
xlabel('Discretization of M_t (mm)','fontsize',15);
set(gca,'fontsize',15,'ytick',[]);
saveas(h1,'figures/ConstCrackDiscrete','fig');
print(h1,'-depsc','figures/ConstCrackDiscrete.eps');

%% Plot the Beta Error
h2=figure('visible','off'); 
hold on;
plot(Lat(1:TotalIter),log10(BetaError_myDD),'--sr','linewidth',2);
plot(Lat(1:TotalIter),log10(BetaError_myDDKL),'--*g','linewidth',2);
plot(Lat(1:TotalIter),log10(BetaError_DDNeil),'--om','linewidth',2);
plot(Lat(1:TotalIter),log10(BetaError_Unif),'--<b','linewidth',2);
legend('RPDD','RPDD with KL error','Neil''s DD, no merge','Uniform');
xlabel('Number of intervals for A_t','fontsize',15)
ylabel('log_1_0(\epsilon_\beta)','fontsize',15);
set(gca,'fontsize',15)
grid on 
saveas(h2,'figures/ConstCrackBetaError','fig');
print(h2,'-depsc','figures/ConstCrackBetaError.eps');

%% Plot the KL distance Error
h3=figure('visible','off'); hold on;
plot(Lat(1:TotalIter),log10(KLError_myDD),'--sr','linewidth',2);
plot(Lat(1:TotalIter),log10(KLError_myDDKL),'--*g','linewidth',2);
plot(Lat(1:TotalIter),log10(KLError_DDNeil),'--om','linewidth',2);
plot(Lat(1:TotalIter),log10(KLError_Unif),'--<b','linewidth',2);

legend('RPDD','RPDD with KL error','Neil''s DD, no merge','Uniform');
xlabel('Number of intervals for A_t','fontsize',15)
ylabel('log_1_0(\epsilon_K_L)','fontsize',15);
set(gca,'fontsize',15)
grid on 
saveas(h3,'figures/ConstCrackKLError','fig');
print(h3,'-depsc','figures/ConstCrackKLError.eps');


%% Plot the mean and standard deviation
h4=figure('visible','off'); hold on;
plot(1:T,mu_da_BLR(2:end),'-k','linewidth',2);
plot(1:T,mu_da_myDD,'--r','linewidth',2);
plot(1:T,mu_da_myDDKL,'--g','linewidth',2);
plot(1:T,mu_da_DDNeil,'--m','linewidth',2);
plot(1:T,mu_da_Unif,'--b','linewidth',2);
axis([0 20 0.9 1.2])
legend('BLR','RPDD','RPDD with KL error','Neil''s DD, no merge','Uniform');
xlabel('Time step','fontsize',17);
ylabel('Mean of dA','fontsize',17);
set(gca,'fontsize',17);
saveas(h4,'figures/ConstCrack_daMean','fig');
print(h4,'-depsc','figures/ConstCrack_daMean.eps');

h5=figure('visible','off'); hold on;
plot(1:T,sqrt(diag(sigma2_da_BLR(2:end,2:end))),'-k','linewidth',2);
plot(1:T,std_da_myDD,'--r','linewidth',2);
plot(1:T,std_da_myDDKL,'--g','linewidth',2);
plot(1:T,std_da_DDNeil,'--m','linewidth',2);
plot(1:T,std_da_Unif,'--b','linewidth',2);
axis([0 20 0.285 0.311])
legend('BLR','RPDD','RPDD with KL error','Neil''s DD, no merge','Uniform','Location','southeast');
xlabel('Time step','fontsize',17);
ylabel('Standard deviation of dA','fontsize',17);
set(gca,'fontsize',17);
saveas(h5,'figures/ConstCrack_daStd','fig');
print(h5,'-depsc','figures/ConstCrack_daStd.eps');


h6=figure ('visible','off'); hold on;
plot(1:T,mu_at_BLR,'-k','linewidth',2);
plot(1:T,mu_at_myDD,'--r','linewidth',2);
plot(1:T,mu_at_myDDKL,'--g','linewidth',2);
plot(1:T,mu_at_DDNeil,'--m','linewidth',2);
plot(1:T,mu_at_Unif,'--b','linewidth',2);

legend('BLR','RPDD','RPDD with KL error','Neil''s DD, no merge','Uniform','Location','southeast');
xlabel('Time step','fontsize',17);
ylabel('Mean of A_t','fontsize',17);
set(gca,'fontsize',17);
saveas(h6,'figures/ConstCrack_atMean','fig');
print(h6,'-depsc','figures/ConstCrack_atMean.eps');


h7=figure('visible','off'); hold on;
plot(1:T,sqrt(sigma2_at_BLR),'-k','linewidth',2);
plot(1:T,std_at_myDD,'--r','linewidth',2);
plot(1:T,std_at_myDDKL,'--g','linewidth',2);
plot(1:T,std_at_DDNeil,'--m','linewidth',2);
plot(1:T,std_at_Unif,'--b','linewidth',2);

legend('BLR','RPDD','RPDD with KL error','Neil''s DD, no merge','Uniform');
xlabel('Time step','fontsize',17);
ylabel('Standard deviation of A_t','fontsize',17);
set(gca,'fontsize',17);
saveas(h7,'figures/ConstCrack_atStd','fig');
print(h7,'-depsc','figures/ConstCrack_atStd.eps');


%% Plot the failure probability
h8=figure('visible','off'); hold on;
plot(1:T,Beta_BLR,'-k','linewidth',2);
plot(1:T,Beta_myDD,'--r','linewidth',2);
plot(1:T,Beta_myDDKL,'--g','linewidth',2);
plot(1:T,Beta_DDNeil,'--m','linewidth',2);
plot(1:T,Beta_Unif,'--b','linewidth',2);

legend('BLR','RPDD','RPDD with KL error','Neil''s DD, no merge','Uniform');
xlabel('Time step','fontsize',15)
ylabel('\beta','fontsize',15)
set(gca,'fontsize',15);
grid on
saveas(h8,'figures/ConstCrackBeta','fig');
print(h8,'-depsc','figures/ConstCrackBeta.eps');

















