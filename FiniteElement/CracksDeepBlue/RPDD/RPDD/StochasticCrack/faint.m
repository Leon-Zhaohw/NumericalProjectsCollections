function at=faint(atm,Yt,C,m,DeltaS,DeltaN)
F=GeometricFun(atm);
if (1-m/2)*C*Yt.*(F*DeltaS)^m*pi^(m/2)*DeltaN+atm.^(1-m/2)<0
    %?????????
    at=1000;
    return
end
at=((1-m/2)*C*Yt.*(F*DeltaS)^m*pi^(m/2)*DeltaN+atm.^(1-m/2)).^(1/(1-m/2));
