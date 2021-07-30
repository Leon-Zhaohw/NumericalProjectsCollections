function result = ComputeTensorEntry( i1, i2, j1,j2, k1, k2)
%Compute the tensor entries for one Neumann entries
%   Detailed explanation goes here

fun1 = @(x)cos(i1*x).*sin(j1*x).*sin(k1*x);
val1 = integral(fun1,0,pi);

fun2 = @(x)cos(j2*x).*sin(i2*x).*sin(k2*x);
val2 = integral(fun2,0,pi);

fun3 = @(x)cos(j1*x).*sin(i1*x).*sin(k1*x);
val3 = integral(fun3,0,pi);

fun4 = @(x)cos(i2*x).*sin(j2*x).*sin(k2*x); 
val4 = integral(fun4,0,pi);

root_inv_lambda = 1.0 / sqrt(i1*i1 + i2*i2);

result = (i1*j2*val1*val2 - i2*j1*val3*val4) * root_inv_lambda * 4 / pi / pi;

end