function [X_mean, X_Std]=MeanStd(Xmarg, X_R)

X_Rmid=[(X_R(2:end-1)+X_R(1:end-2))/2 X_R(end-1)];

X_mean=sum(X_Rmid.*Xmarg');

X_Std= sqrt(sum(X_Rmid.^2.*Xmarg')-X_mean^2);
