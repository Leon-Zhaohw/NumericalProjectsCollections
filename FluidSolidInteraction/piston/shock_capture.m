% -------------------------------
% filename shock_capture.m
%
% if ifct is set to zero we compute a solution of order 1
% if ifct is set to 2 we compute a solution of order 2
% if ifct is set to 1 we combine both low and high order solutions
%

function[du]=shock_capture(ifct,vcor,conec,vmgn,vmgnp1,vres,vsol,xlumpm,number);
%
% we recall that :
% nnt  number of nodes (i.e. discretization points)
% ndln number of dof per node
% nelt number of finite element for C-model
% nnel number of node per finite element
%
[nnt  ndln] = size(vsol);
[nelt nnel] = size(conec);
%
% Numerical diffusive coefficient for first order solution.
%
cd   = 0.4;
%
% intermediate tables used to locate numerical oscillation
%
pmax = zeros(nnt,1);
pmin = zeros(nnt,1);
rmax = zeros(nnt,1);
rmin = zeros(nnt,1);
%
% Loop over the nodal dof (volumic mass, momentum, total energy per volum unit)
%
for j=1:ndln
    vresj = vres(:,j);
    vsolj = vsol(:,j);
    %
    % Correction of the residual due to the moving mesh
    %
    vresj = vresj+(vmgn-vmgnp1)*vsolj;
    %
    % Low order solution (first order)
    %
    vdul = (vresj+cd.*vmgn*vsolj -cd.*xlumpm.*vsolj)./xlumpm;
    %
    % High order solution (second order)
    %
    vduh = vmgnp1\vresj;
    %
    vdul3(:,j) = vdul;
    vduh3(:,j) = vduh;
    %
    % Shock capturing technique
    %
    if(ifct~=0)
        uimax = -1.0E10.*ones(nnt,1);
        uimin = +1.0E10.*ones(nnt,1);
        %
        % loop over the finite elements
        %
        for ie=1:nelt
            kloce = conec(ie,1:nnel);
            vcore = vcor(kloce,:); x = vcore(2 ,:) - vcore(1,:); xl=sqrt(x * x' );

            uimax(kloce)  = max(uimax(kloce),...
                max(max(vsolj(kloce),(vsolj(kloce)+vdul(kloce)))));
            uimin(kloce)  = min(uimin(kloce), ...
                min(min(vsolj(kloce),(vsolj(kloce)+vdul(kloce)))));
            aec(ie,[1 2]) = ((xl./6.*[1 -1;-1 1]*(cd.*vsolj(kloce)+vduh(kloce)))...
                ./xlumpm(kloce))';
            pmax(kloce) = pmax(kloce)+max(zeros(2,1),aec(ie,[1 2])');
            pmin(kloce) = pmin(kloce)+min(zeros(2,1),aec(ie,[1 2])');
        end

        qmax = uimax-(vsolj+vdul);
        qmin = uimin-(vsolj+vdul);
        %
        % loop over the nodes
        %
        for in=1:nnt
            if (pmin(in)==0.) | (pmax(in)==0.)
                rmax(in)=0. ; rmin(in)=0.;
            else
                rmax(in)=min(1.,qmax(in)./pmax(in));
                rmin(in)=min(1.,qmin(in)./pmin(in));
            end
        end
        %
        % loop over the finite elements
        %
        for ie=1:nelt
            for in=1:2
                if (aec(ie,in))>=0.
                    cel(ie,in)=min(rmax(conec(ie,1)),rmax(conec(ie,2)));
                else
                    cel(ie,in)=min(rmin(conec(ie,1)),rmin(conec(ie,2)));
                end
            end
        end
        %
        if j==1 cel1=cel;  end;
        if j==2 cel2=cel;  end;
        if j==3 cel3=cel;  end;
    end
end

var_cel=zeros(nnt,1);
%
% We compute the slope limiter (switch between first and second order)
%
switch ifct
    case 0
        var_cel=zeros(nnt,1);
    case 1
        for ie=1:nelt
            celi=min(cel1(ie,:),min(cel2(ie,:),cel3(ie,:)));
            var_cel(conec(ie,:))=var_cel(conec(ie,:))+celi'./number(conec(ie,:));
        end
    case 2
        var_cel=ones(nnt,1);
end
%
% Loop over the nodal dof (volumic mass, momentum, total energy per volum unit)
%
% we update the increment step
%
for j=1:ndln
    du(:,j) = vdul3(:,j)+var_cel.*(vduh3(:,j)-vdul3(:,j));
end
