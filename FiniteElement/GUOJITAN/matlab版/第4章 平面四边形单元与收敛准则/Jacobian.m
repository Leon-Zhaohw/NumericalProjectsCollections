function [N_cxi, N_eta, JA, dt_JA, inv_JA] = Jacobian(X0,Y0 )
% 根据四个结点坐标， 推算雅克比矩阵、雅克比行列式、雅克比逆阵
%输入参数：4个结点坐标值X0(1:4,1),Y0(1:4,1)列向量
%输出量：带有符号的表达式—形函数对局部坐标导数、雅克比矩阵、雅克比行列式、雅克比逆阵
syms  cxi  eta             %定义符号参量：cxi——  ， eta——   ；
%  计算形函数对局部坐标的导数
N_cxi = [-(1-eta),  (1-eta),  (1+eta),  -(1+eta)]/4;     %形函数对无量纲局部坐标偏导数
N_eta =[-(1-cxi),  -(1+cxi),  (1+cxi),  (1-cxi)]/4;
%  计算雅克比矩阵
JA= [ N_cxi*X0 ,  N_cxi*Y0 ;
    N_eta*X0,   N_eta*Y0]
%  计算雅克比行列式
dt_JA=det(JA)
%  计算雅克比逆阵
inv_JA=inv(JA)
return


