function  Yt_cpd=fYint(alphaY,Ytm,Yt,Y_mean,Y_std,Y_mu,Y_sigma,DeltaN) 

Z_sigma=Y_sigma;
rouz=(log(Y_std^2/Y_mean^2*exp(-DeltaN./alphaY)+1))/Y_sigma^2;
Ztc_sigma=Z_sigma.*sqrt(1-rouz.^2);
Ztm=log(Ytm)-Y_mu;
Ztc_mu=rouz.*Ztm;
Zt=log(Yt)-Y_mu;
Yt_cpd=0.5*erf((Zt-Ztc_mu)./(sqrt(2)*Ztc_sigma))+0.5;
