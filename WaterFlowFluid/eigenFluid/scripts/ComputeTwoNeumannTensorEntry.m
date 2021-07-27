function result = ComputeTwoNeumannTensorEntry( i1, i2, j1,j2, k1, k2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fun1 = @(x)sin(i1*x).*cos(j1*x).*cos(k1*x);
val1 = integral(fun1,0,pi);

fun2 = @(x)sin(i2*x).*cos(j2*x).*sin(k2*x);
val2 = integral(fun2,0,pi);

fun3 = @(x)cos(i1*x).*sin(j1*x).*cos(k1*x);
val3 = integral(fun3,0,pi);

fun4 = @(x)cos(i2*x).*sin(j2*x).*sin(k2*x);
val4 = integral(fun4,0,pi);

root_inv_lambda = 1.0 / sqrt(i1*i1 + i2*i2);

if (i2 == 0)
    root_inv_lambda = 0;
end

if (i1 == 0 && i2 > 0)
    root_inv_lambda = 1.0 / sqrt(2) / i2; 
end

n_factor = 4 / pi / pi;

%both zero
if (k1 == 0 && j1 == 0)
    n_factor = n_factor * 0.5;
elseif (k1 == 0 || j1 == 0)  %one zero
    n_factor = n_factor / sqrt(2);
end

result = (i1*j2*val1*val2 - i2*j1*val3*val4) * root_inv_lambda * n_factor;

end