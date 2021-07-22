function [x,y]=XY_CxiEta(X0,Y0)
% 根据四个结点坐标， 推算四点等参元的几何模式
% 输入参数：4个结点坐标值X0(1:4,1),Y0(1:4,1)列向量
%输出量：带有符号的表达式——直角坐标系下的坐标x、y的表达式
syms  cxi  eta              %符号：cxi——  ， eta——   ；
N(1)=(1-cxi)*(1-eta)/4;      %形函数N(1,1:4)，行向量；
N(2)=(1+cxi)*(1-eta)/4;
N(3)=(1+cxi)*(1+eta)/4;
N(4)=(1-cxi)*(1+eta)/4;
x=N*X0;      y=N*Y0;             %按式（4.31）计算
x=expand(x);                 %将x、y用普通多项式表示
y=expand(y);
return

