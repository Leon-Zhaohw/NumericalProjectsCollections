result = zeros(1,100);
for i = 1:100
    fun1 = @(x)sin((i-0.5)*x).*sin(i*x);
    val1 = integral(fun1,0,pi);
    result(i) = val1;
end

plot(result);