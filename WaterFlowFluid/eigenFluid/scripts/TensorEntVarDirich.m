function result = TensorEntVarDirich( i1, i2, j1,j2, k1, k2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
lambda_i = sqrt(i1*i1 + i2*i2);
inv_lambda_j = 1.0 / sqrt(j1*j1 + j2*j2);
inv_lambda_k = 1.0 / sqrt(k1*k1 + k2*k2);
c = lambda_i * inv_lambda_j * inv_lambda_k * 8 / (pi^3);

fun1 = @(x)sin(i1*x).*cos(j1*x).*sin(k1*x);
val1 = integral(fun1,0,pi);
fun2 = @(x)sin(i2*x).*sin(j2*x).*cos(k2*x);
val2 = integral(fun2,0,pi);
fun3 = @(x)sin(i1*x).*sin(j1*x).*cos(k1*x);
val3 = integral(fun3,0,pi);
fun4 = @(x)sin(i2*x).*cos(j2*x).*sin(k2*x);
val4 = integral(fun4,0,pi);

result = (-j1*k2*val1*val2 + j2*k1*val3*val4)*c;

end

