

function [Q1,Q2]=Equivalent_Nodal_Force_Surface(i1,i2,dd,QD,kd)
%计算线性分布面力，边界上二结点的面力等效载荷
%输入量：i1,i2边界二个结点号；dd,QD有面力的结点数及面力值数据；
% kd面力类型代号：1——法向面力，2——切向面力，3——面力分量表示的面力
%返回值：二个结点的面力等效载荷大小，2个向量（2*1）
global  t0  XY                 %t0板厚，结点坐标
Q1=zeros(2,1); Q2=zeros(2,1);
if kd<3
    qi=0; qj=0;                    %法向、切向面力将两端面力清零
    for t=1:dd                        %dd——面力的结点个数
        if QD(t,1)==i1
            qi= QD(t,3);                 %提取对应线段 1#结点的面力值
        elseif QD(t,1)==i2
            qj= QD(t,3);                %提取对应线段 1#结点的面力值
        end
    end
else
    qi=zeros(2,1); qj=zeros(2,1);         %坐标分量两端面力二个方向清零
    for t=1:dd
        if QD(t,1)==i1
            qi(1:2,1)= QD(t,3:4);                % 1#结点面力：二个方向分量
        elseif QD(t,1)==i2
            qj(1:2,1)= QD(t,3:4);               % 1#结点面力：二个方向分量
        end
    end
end
p1= t0* (2*qi+qj)/6;
p2= t0* (qi+2*qj)/6;
ax=XY(i2,1)-XY(i1,1); ay=XY(i2,2)-XY(i1,2);
switch kd
    case 1                            %法向面力
        BL=[-ay; ax];
    case 2                            %切向面力
        BL=[ax; ay];
    case 3                            %斜向面力
        BL=sqrt(ax*ax+ay*ay);
end
Q1= p1*BL;        % 1#结点的等效载荷
Q2= p2*BL;        % 2#结点的等效载荷
return



