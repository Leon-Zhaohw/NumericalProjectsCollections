function Epsilon_X=KLError(X_R,Xmarg)
LX=length(Xmarg);
Epsilon_X=zeros(1,LX);
for i=1:LX
    if i==1
        fmin=1E-10;
        fmax=max([Xmarg(i) Xmarg(i+1)]);
    elseif i==LX
        fmax=max([Xmarg(i-1) Xmarg(i)]);
        fmin=1E-10;
    else
        fmax=max([Xmarg(i-1) Xmarg(i+1) Xmarg(i)]);
        fmin=min([Xmarg(i-1) Xmarg(i+1) Xmarg(i)]);
    end
    fmid=Xmarg(i);
    fbar=(fmax+fmin+fmid)/3;
    Epsilon_X(i)=((fmax-fbar)/(fmax-fmin)*fmin*log(fmin/fbar)+(fbar-fmin)/(fmax-fmin)*fmax*log(fmax/fbar))*abs(X_R(i+1)-X_R(i));
    if isnan(Epsilon_X(i)) || isinf(Epsilon_X(i))
        Epsilon_X(i)=0;
    end
end
