function F=GeometricFun(atm)
%Edge crack case
PlateWidth=100;
alpha=atm/PlateWidth;
F=0.265*(1-alpha).^4+(0.857+0.265*alpha)./(1-alpha).^1.5;
if imag(F)~=0
    disp('Warning: complec geometric function')
end