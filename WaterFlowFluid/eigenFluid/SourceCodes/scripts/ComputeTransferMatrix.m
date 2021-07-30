T = zeros(20,20);
for k = 1:20
    for i = 1:20
        T(k,i) = 0.5 * ( sin((k-i+0.5)*pi)/(k-i+0.5) - sin((k+i-0.5)*pi)/(k+i-0.5)) * 2 / pi;
    end
end


[Q,R] = qr(T);
diag(transpose(Q) * Q);
figure;
image(T*255);
figure;
image(-Q*255);