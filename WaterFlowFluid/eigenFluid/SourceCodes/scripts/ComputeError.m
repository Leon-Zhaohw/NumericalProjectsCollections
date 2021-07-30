n = 40;
C = zeros(n,n);
for k  = 1:n
    for i = 1:n
        C(k,i) = 0.5 * (sin((k-i+0.5)*pi)/(k-i+0.5) - sin((k+i-0.5)*pi)/(k+i-0.5));
    end
end

Err = zeros(n,1);
x = 0:pi/10000:pi;

for k = 1:n
    Approx = zeros(1,10001);
    for i = 1:n
        Approx = Approx + C(k,i)*sin((i-0.5)*x)*2 / pi;
    end
    EV = sin(k*x) - Approx;
    Err(k) = sum(EV.*EV)/10000;
end

figure
plot(Err)
%xlabel('k')
%ylabel('Error magnitude')