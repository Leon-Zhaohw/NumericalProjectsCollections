function [mu,sigma]=LNormToNorm(mean,std)
%Transform the mean and std in lognormal to mu and sigam
mu=log(mean)-1/2*log(1+(std/mean)^2);
sigma=(log(1+(std/mean)^2))^0.5;