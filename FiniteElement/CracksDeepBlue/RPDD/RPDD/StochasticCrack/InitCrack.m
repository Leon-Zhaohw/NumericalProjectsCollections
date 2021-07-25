function [a1_R,a1cdf]=InitCrack(at_R,atmarg)


a1cdf=zeros(1,length(atmarg));
a1cdf(1)=atmarg(1);
for i=2:length(a1cdf)
    a1cdf(i)=a1cdf(i-1)+atmarg(i,1);
end
a1_R=at_R(1:end-1);