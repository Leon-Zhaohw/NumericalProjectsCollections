

function GI = Integral_Gauss(f,a,b,n)
%高斯积分—完成1~5个高斯点的数值积分；若超出该范围予以提示并结束。
% f—被积函数，表达式。
%输入量：a、b积分区间[a,b]上限、下限；n积分点个数。
if(n>5 || n<1)
    fprint('您输入的高斯积分点数n=','%4i','未在[1,5]范围内，请修改后再运行。',n)
    return
end
ta = (b-a)/2;
tb = (a+b)/2;
switch n
    case 1,
        ci=[0];
        wi=[2];
    case 2,
        ci=[-0.5773503, 0.5773503];
        wi=[1,1];
    case 3,
        ci=[-0.7745967,0 ,0.7745967];
        wi=[0.55555556, 0.88888889,0.55555556];
    case 4,
        ci=[-0.8611363, -0.3398810, 0.3398810 , 0.8611363];
        wi=[0.3478548, 0.6521452, 0.6521452, 0.3478548];
    case 5,
        ci=[-0.9061793, -0.5384693, 0 , 0.5384693, 0.9061793];
        wi=[0.2369269, 0.4786287, 0.5688889, 0.4786287, 0.2369269];
end
GI=0;
for j=1:n
    GI=GI+ta*(wi(j)*subs(sym(f),findsym(sym(f)),(ta*ci(j)+tb)));
end
return
