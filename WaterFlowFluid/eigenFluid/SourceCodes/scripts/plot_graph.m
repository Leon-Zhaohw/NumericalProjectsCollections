function plot_graph( C,k )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Approx = zeros(1,10001);
x = 0:pi/10000:pi;
Colnum = size(C,2);
for i = 1:20
    Approx = Approx + C(k,i)*sin((i-0.5)*x) * 2/pi;
end
y = sin(k*x);
Erf = y -Approx;
yg = sin(1.5*x);
%figure;
%plot(Erf);

yy = zeros(1,10001);
for i=21:21
    yy = yy + C(k,i)*sin((i-0.5)*x) * 2/pi;
end

plot(x, yy, x , Erf, '--');
sum(yg.*yy)
%figure;
%plot(Dotval);

end

