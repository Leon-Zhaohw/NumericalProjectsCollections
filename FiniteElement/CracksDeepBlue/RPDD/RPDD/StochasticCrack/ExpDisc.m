function CPD=ExpDisc(x,mean)
Fx=expcdf(x,mean);
Fx1=Fx(2:end);
Fx2=Fx(1:end-1);
CPD=Fx1-Fx2;
