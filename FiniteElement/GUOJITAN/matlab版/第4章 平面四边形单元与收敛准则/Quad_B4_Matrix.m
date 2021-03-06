
function [B, dt] = Quad_B4_Matrix( ie, cxi, eta )
%  计算四点等参元几何矩阵B在某点的值
%  输入参数：ie ——单元号，cxi,eta ——位置点局部坐标
%  返回量：B——对直角坐标系的几何矩阵B
%          dt——det(JA)雅克比行列式
global  EL  XY             %全局变量，EL——单元信息， X，Y——结点坐标
X0(1:4,1)=XY(EL(ie,1:4),1);   %将单元ie的4个结点x坐标值，赋予列向量X0(4)
Y0(1:4,1)=XY(EL(ie,1:4),2);    % y坐标值，赋予列向量Y0(4)
N_cxi = [-(1-eta),  (1-eta),  (1+eta),  -(1+eta)]/4;  %形函数对无量纲局部坐标偏导数
N_eta =[-(1-cxi),  -(1+cxi),  (1+cxi),  (1-cxi)]/4;
%  计算雅克比矩阵
JA =  [ N_cxi*X0,  N_cxi*Y0
    N_eta*X0,   N_eta*Y0 ];
%  计算雅克比行列式
dt=det(JA);
%  计算雅克比逆阵
in_JA=inv(JA);
%  计算形函数对直角坐标系的偏导数
A = in_JA*[N_cxi; N_eta] ;     %中括号内为2×4，第一行为N_cxi;第二行为 N_eta
N_x = A(1,:) ;
N_y = A(2,:) ;
% 组装几何矩阵
B = zeros( 3, 8 ) ;
for i=1:1:4
    B(1:3,(2*i-1):2*i) = [ N_x(i) ,     0 ;      %B中1~3行，(2*i-1)~2*i列的子阵
        0,   N_y(i);
        N_y(i),    N_x(i) ] ;
end
return

















