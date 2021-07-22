

function [B, A] =  Plane_B3_Matrix( ie)
%  计算单元几何矩阵B、面积A
%  输入参数:ie ----  单元号,
%  返回量：几何矩阵B和三角形单元面积A
global  XY  EL       %全局变量：结点坐标、单元信息
i=EL(ie,1);j=EL(ie,2);m=EL(ie,3);       %单元的3个结点
bi = XY(j,2) - XY(m,2);
bj = XY(m,2) - XY(i,2);
bm = XY(i,2) - XY(j,2);
ci = XY(m,1) - XY(j,1);
cj = XY(i,1) - XY(m,1);
cm = XY(j,1) - XY(i,1);
A = (bj * cm - bm * cj )/2;            %单元面积
B = [bi  0  bj  0  bm  0 ;               %几何矩阵
    0   ci  0  cj  0  cm ;
    ci  bi  cj  bj  cm  bm]/(2*A);
return

