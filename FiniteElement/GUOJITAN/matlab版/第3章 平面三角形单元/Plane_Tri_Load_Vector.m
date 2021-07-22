

function ZQ= Plane_Tri_Load_Vector
%计算平面问题各种外力的等效结点载荷，且组装到总载荷向量；
%外力的形式包括：竖向重力、不等的体力、结点集中力、线性分布的法向、切向、斜向面力等
%返回量：总载荷向量ZQ
global  t0  lou  nd  ne  nj   nt  XY  EL  TL QJ
%全局变量： t0,lou --厚度、密度； nd ne ng nj nt--结点数,单元数,约束总数,集中力个数,体力类型数
%  XY、EL、TL、QJ--数组： 节点坐标、单元信息、各类体力数值、集中力信息
global  mx md QMD QMX
%   三类面力的：有面力线段数、有结点面力个数、各点面力值、各段的两端结点及面力类型
ZQ = zeros(2*nd,1 ) ;      % 结构的总载荷列阵

%step1：计算体力的等效结点载荷，并组装到总的载荷向量中
if lou~=0  |    nt>0                %密度不为0或存在体力时，计算重力或体力的等效结点载荷
    for ie=1:1:ne
        i=EL(ie,1);j=EL(ie,2);m=EL(ie,3);
        xi=XY(i,1);xj=XY(j,1);xm=XY(m,1);
        yi=XY(i,2);yj=XY(j,2);ym= XY(m,2);
        A = (xi*(yj-ym) + xj*(ym-yi) + xm*(yi-yj))/2;
        if lou~=0                               %计算重力
            qy = - A*t0*lou/3;
            for s=1:1:3
                i0= 2*EL(ie,s);
                ZQ(i0,1) = ZQ(i0,1) + qy;
            end
        end
        if  nt>0                                %存在体力
            it0=EL(ie,4);
            Qx=TL(it0,:)*A*t0/3;
            for s=1:1:3
                i0= 2*EL(ie,s);
                ZQ(i0-1,1) = ZQ(i0-1,1) +  Qx(1);
                ZQ(i0,1) = ZQ(i0,1) +  Qx(2);
            end
        end
    end
end

% step2：将结点集中力组装到总的载荷向量中
if nj>0
    for s=1:1:nj
        i0=2*QJ(s,1);
        ZQ(i0-1,1) = ZQ(i0-1,1) + QJ(s,2);
        ZQ(i0,1) = ZQ(i0,1) + QJ(s,3);
    end
end
% step3：计算线性分布压力的等效结点载荷，并组装到总的载荷向量中
if  mx>0
    mx1=0;mx2=0;mx3=0;
    index = (QMX(:,3) ==1);         %将面力分类，有法向面力作用的边界
    if  ~isnan(index)
        QMX1=QMX(index,:);
        [mx1,m0]=size(QMX1);              %法向面力作用边界的数量
        row_index =( QMD(:,2) ==1);          %法向面力的结点
        if ~isnan(row_index)
            QMD1=QMD(row_index,:);
            [md1,m0]=size(QMD1);               %有法向面力值的结点数量
        else
            disp('面力数据存在错误，缺少法向面力的结点值')
        end
    end
    index = (QMX(:,3) ==2);              %有切向面力作用的边界
    if  ~isnan(index)
        QMX2=QMX(index,:);
        [mx2,m0]=size(QMX2);               %切向面力作用边界的数量
        row_index =( QMD(:,2) ==2);             %切向面力
        if ~isnan(row_index)
            QMD2=QMD(row_index,:);
            [md2,m0]=size(QMD2);               %有切向面力值的结点数量
        else
            disp('面力数据存在错误，缺少切向面力的结点值')
        end
    end
    
    index = (QMX(:,3) ==3);               %存在以坐标分量表示的斜向面力
    if  ~isnan(index)
        QMX3=QMX(index,:);
        [mx3,m0]=size(QMX3);               %存在以坐标分量表示斜向面力的数量
        row_index =( QMD(:,2) ==3);           %以坐标分量表示的斜向面力
        if ~isnan(row_index)
            QMD3=QMD(row_index,:);
            [md3,m0]=size(QMD3);          %有斜向（体力分量）面力值的结点数量
        else
            disp('面力数据存在错误，缺少以面力分量表示的结点值')
        end
    end
    
    if  mx1>0
        for s=1:1:mx1
            i1=QMX1(s,1);i2= QMX1(s,2);
            [Q1,Q2]=Equivalent_Nodal_Force_Surface(i1,i2,md1,QMD1,1);
            ZQ([2*i1-1,2*i1],1) =  ZQ([2*i1-1,2*i1],1) + Q1;
            ZQ([2*i2-1,2*i2],1) =  ZQ([2*i2-1,2*i2],1) + Q2;
        end
    end
    
    if  mx2>0
        for s=1:1:mx2
            i1=QMX2(s,1);i2= QMX2(s,2);
            [Q1,Q2]=Equivalent_Nodal_Force_Surface(i1,i2,md2,QMD2,2);
            ZQ([2*i1-1,2*i1],1) =  ZQ([2*i1-1,2*i1],1) + Q1;
            ZQ([2*i2-1,2*i2],1) =  ZQ([2*i2-1,2*i2],1) + Q2;
        end
    end
    
    if  mx3>0                   % nx3----有斜向线性分布表面力的线段个数
        for s=1:1:mx3
            i1=QMX3(s,1);i2= QMX3(s,2);
            [Q1,Q2]=Equivalent_Nodal_Force_Surface(i1,i2,md3,QMD3,3);
            ZQ([2*i1-1,2*i1],1) =  ZQ([2*i1-1,2*i1],1) + Q1;
            ZQ([2*i2-1,2*i2],1) =  ZQ([2*i2-1,2*i2],1) + Q2;
        end
    end
end
return


