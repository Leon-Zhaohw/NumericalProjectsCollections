%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test LGL collocation method for the second-order BVP:
% -u''(x)+u(x)=f(x),
% f=1-x-x^2/2, -1<x<0, f=1-x, 0<=x<1
% with exact solution
% u=cosh(x+1)-x^2/2-x, -1<x<0; u=cosh(x+1)-cosh(x)-x+1, 0<=x<1.
% u\in C^3([0,1]), f\in C^1([0,1])
% Use: legslb()--LGL nodes & weights; legslbdiff()---1st-order Diff. matrix
% Create on June 22, 2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clf; format long e
nv=[8 16 32 64 96 128 256 512 800 1024]'; % Number of collocation (LGL) points   
cmerr=[]; % Store the maximum pointwise error for N
for l=1:length(nv)
    %% Initialization: LGL nodes & weights; Legendre transform matrices 
    N=nv(l); % N+1 collocation points 
    [x,w]=legslb(N+1);   % Compute LGL points (in ascending order) x, and weights w
    %% Exact solution and source term
    ue=zeros(N+1,1); % Initialize the exact solution and its derivatives     
    f=zeros(N+1,1); % Initialize the RHS f(x)
    k=find(x<0,1,'last'); % find the last index such that x<0
    x1=x(1:k); x2=x(k+1:end); % x1<0; x2>0   
    ue(1:k)=cosh(x1+1)-x1.^2/2-x1; ue(k+1:N+1)=cosh(x2+1)-cosh(x2)-x2+1; % exact solution
    f(1:k)=1-x1.^2/2-x1; f(k+1:N+1)=1-x2; 
    cminus=ue(1); cplus=ue(N+1);  % Boundary conditions taken from exact solution 
    %% Collocation Method
    D=legslbdiff(N+1,x); % Compute 1st-order Differentiation Matrix
    D2=D*D; D2in=D2(2:N,2:N); % 2nd-order Diff matrix and its interior part
    A=-D2in+eye(N-1); % Coefficient matrix of 
    uB=-cminus*D2(2:N,1)-cplus*D2(2:N,N+1); % RHS: induced by the B.C.   
    fv=f(2:N)-uB; % RHS vector
    uN=A\fv; % Numerical solution 
    cmerr=[cmerr, max(abs(ue(2:N)-uN))]; % Max-pointwise error
end

%% plot exact and numerical solution at N=128
 figure(1)
 plot(x,ue,'r',x,[cminus;uN;cplus],'b-.','Markersize',8.0,'Linewidth',2.0);
 set(gca,'xtick',-1:0.5:1,'ytick',1.3:0.1:2.3,'fontsize',16)
 title('Exact Solution Vs Numerical Solution')
 xlabel('x','fontsize',14); ylabel('u(x) & u_N(x)','fontsize',14)
 legend('u(x)','u_N(x)','location','NorthWest')
 grid on

%% plot numerical errors
 figure(2)
 plot(log10(nv),log10(cmerr),'-.o','Markersize',8.0,'Linewidth',2.0);
 hold on
 plot(log10(nv),log10(nv.^(-4)),'-.','Markersize',10.0,'Linewidth',2.0);
 xlabel('log_{10}(N)','fontsize',14); ylabel('log_{10}(error)','fontsize',14)
 axis([0.5 3.5 -13 -2])
 set(gca,'xtick',0.5:0.5:3.5,'ytick',-14:2:0,'fontsize',16)
 legend('max|u(x_j)-u_N(x_j)|','N^{-4}','location','NorthEast')
 title('Maximum Pointwise Error')
 grid on
