function [alphaY_CPD  Y1_CPD Ytc_CPD a1_CPD,at_CPD,M_CPD]=CPDgen(alphaY_R,Y_R,at_R,M_R,alphaY_mean,alphaY_std,Y_mean,Y_std,a0_mean,M_sigma,DeltaN,DeltaS,m,C,ac)


[alphaY_mu alphaY_sigma]=LNormToNorm(alphaY_mean,alphaY_std);
alphaY_CPD=NormalDisc(log(alphaY_R),alphaY_mu,alphaY_sigma);
LalphaY=length(alphaY_R)-1;


%Y(n) is stationary lognormal process
[Y_mu Y_sigma]=LNormToNorm(Y_mean,Y_std);
Y_CPD=NormalDisc(log(Y_R),Y_mu,Y_sigma);
LY=length(Y_R)-1;
Y1_CPD=zeros(LalphaY,LY);
for i=1:LalphaY
    for j=1:LY
        Y1_CPD(i,j)=Y_CPD(j);
    end
end

Y_Rinf=8;
lambda=10;

Ytc_CPD=Ytc_CPDgen(alphaY_R,Y_R,Y_Rinf,Y_mean,Y_std,Y_mu,Y_sigma,DeltaN);

a0_CPD=ExpDisc(at_R,a0_mean);

%small lambda should be a large Y_Rinf value
Y_Rinf=10;   % This limit of Y_Rinf is diffrent due to the sigularity point in the range
[at_CPD,a1_CPD]=a_CPDgen(a0_CPD,at_R,Y_R,C,m,DeltaS,DeltaN);

lambda=1;
a_Rinf=500;
M_CPD=M_CPDgen(at_R,M_R,M_sigma,lambda,ac,a_Rinf);
disp('Finish generating CPD')

