function [FailProb_BLR mu_at_BLR sigma2_at_BLR mu_da_BLR sigma2_da_BLR]=ConstCrack_BLR(Evidence,T,mu_a0,sigma_a0,mu_da,sigma_da,sigma_error,atMax)
%The Bayesian linear regression model for exact solution of constant crack
%growth, all variables are assumed to satisfy Gaussian distribution
% disp('******************* Bayesian Linear Regression ********************')
% Evidence=[5 10
%          10 14
%          15 25];
% T=20;
% mu_a0=7;
% sigma_a0=3;
% mu_da=1;
% sigma_da=0.3;
% sigma_error=2;
% atMax=30;

SizeEV=size(Evidence);
X_ev=zeros(T+1,SizeEV(1));
Y_ev=zeros(SizeEV(1),1);
for i=1:SizeEV(1)
      X_ev(1:Evidence(i,1)+1,i)=1;
      Y_ev(i)=Evidence(i,2);
end
sigma2_w=zeros(T+1,T+1);
mu_w=zeros(T+1,1);
for i=1:T+1
    if i==1
        sigma2_w(i,i)=sigma_a0^2;
        mu_w(i)=mu_a0;
    else
        sigma2_w(i,i)=sigma_da^2;
        mu_w(i)=mu_da;
    end
end

A_inv=inv(sigma2_w)+1/sigma_error^2*X_ev*X_ev';
sigma2_da_BLR=inv(A_inv);
mu_da_BLR=sigma2_da_BLR*(inv(sigma2_w)*mu_w+1/sigma_error^2*X_ev*Y_ev);

FailProb_BLR=zeros(1,T);
for t=1:T
    Xpred=zeros(1,T+1);
    Xpred(1:t+1)=1;
    mu_at_BLR(t)=Xpred*mu_da_BLR;
    % if sigma2_at=Xpred*A*Xpred'+sigma_error^2, then it is the covariance
    % matrix for M. if sigma2_Ypred(t)=Xpred*A*Xpred, then it is the
    % covariance matrix for at 
    sigma2_at_BLR(t)=Xpred*sigma2_da_BLR*Xpred';%+sigma_error^2;
    temp=(atMax-mu_at_BLR(t))/sqrt(2*sigma2_at_BLR(t));
    FailProb_BLR(t)=1-0.5*(1+erf(hpf(temp,150)));
end

save('FailProb_BLR','FailProb_BLR','mu_at_BLR','sigma2_at_BLR','mu_da_BLR','sigma2_da_BLR');
