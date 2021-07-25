function [FailProb_KS mu_at_KS sigma2_at_KS]=ConstCrack_KF(Evidence,T,mu_a0,sigma_a0,mu_da,sigma_da,sigma_error,atMax)

% This code was used to generate mean and std of contant crack growth case
% using Kalman Smoother method, as the model is linear gaussian state space
% model, therefore, the exact is available. 

% X(t+1) = F X(t) + B*u(t) noise(Q)
% Y(t) = H X(t) + noise(R)
% T is the time slice
% F, B H is set to be 1
% X(t) is the crack size at time step t
% u(t) is the crack increment at time step t
% noise Q is covariance matrix of da and a0
% noise R is the covariance matrix of Y

disp('******************* Kalman Smooth ********************')
% if the there is missing observation, just set nan at t-th entry of y
y=ones(1,T)*nan;
y(Evidence(:,1))=Evidence(:,2);
ss = 1; % state size
os = 1; % observation size
F = ones(1,T); 
H = ones(1,T);
B=ones(1,T);
Q = zeros(1,T);
R=zeros(1,T);
for i=1:T
    if i==1
        Q(i)=sigma_a0^2+sigma_da^2;
    else
        Q(i)=sigma_da^2;
    end
    R(i)=sigma_error^2;
end
u=ones(1,T)*mu_da;

initx = mu_a0;
initV = sigma_a0^2;
model=1:T;

[mu_at_KF,sigma2_at_KF] = kalman_filter(y, F, H, Q, R, initx, initV,'model',model,'u',u,'B',B);
% xsmooth is the mean value of smoothing results and Vsmooth is the
% covariance of smoothing results
[mu_at_KS, sigma2_at_KS] = kalman_smoother(y, F, H, Q, R, initx, initV,'model',model,'u',u,'B',B);
for t=1:T
    t
    temp=(atMax-mu_at_KS(t))/sqrt(2*sigma2_at_KS(t));
    FailProb_KS(t)=1-0.5*(1+erf(hpf(temp,100)));
end  
FailProb_KS=double(FailProb_KS);

save('FailProb_KS','FailProb_KS','mu_at_KS','sigma2_at_KF');


