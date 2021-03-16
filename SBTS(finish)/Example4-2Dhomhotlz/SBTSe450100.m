%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex Helmholtz equations is 
%                   -\LapLace U + \sigma1 U +i\sigma2 U = f
%its coresspoding complex symmetric linear systems is
%                   ((K + \sigam1*I) +i\sigam2*I)*u = b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SBTS方法的参数选取列表
%m            8          16         32         64          128           256  
%alpha      3.5065     3.4773     3.4730     3.4722       3.4720       3.4720 
%beta       0.5832     0.5840     0.5841     0.5841       0.5841       0.5841
%Umax       1.4380     1.4351     1.4342     1.4340       1.4339       1.4339
%Umin       0.1474     0.0427     0.0114     0.0030       7.5099e-04   1.8774e-04
%IT          20          20          20        20            10          10
%-------------------取两位有效数字---------------------------------------------------------------------------------------
%alpha(IT)  2.13(10)   2.07(11)   2.05(11)   2.04(10)     2.04(10)  
%res        9.8162     2.7354     2.6081     9.7813       9.9150 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%构造矩阵所需要的一些操作及用法
%从一个6x6矩阵中，提取它的第一行元素形成一个6维行向量。A(i,:)行 ,A(:,j)列方法。
%A(i,:)提取A的第i行
%A(:,j)提取A的第j列
clear;                                      %清除上次所生成的内存中的东西
clc;                                        %清屏
sigma1=50;
sigma2=100;
m=input('please input  m=');                %输8入m的值=8,16,32,64,128,256
alpha=input('please input  alpha=');        %输入参数alpha的值，如果输入alpha=0，程序会计算参数alpha的值
h=1/(m+1);                                      %网格步长
n=m*m;                                      %生成矩阵维数
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%三对角矩阵(-1,2-1)
Vm=(1/(h*h))*V;                             
K=kron(speye(m),Vm)+kron(Vm,speye(m));


W=K+sigma1*speye(n);
W=h*h*W;
T=sigma2*speye(n);
T=h*h*T;

p=(W-T)*ones(n,1);
q=(W+T)*ones(n,1);


if(alpha==0&&m~=128&m~=256&&m~=512)
S=W\T;
eigS = eig(full(S));
Umax = max(abs(eigS))
Umin = min(abs(eigS))
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if(alpha==0&&m==128)
S=W\T;
eigS = eig(full(S));
Umax = max(abs(eigS))
Umin = min(abs(eigS))
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if(alpha==0&&m==256)
eigS = eigs(T,W);
Umax = max(abs(eigS))
Umin = 1.8774e-04;
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if(alpha==0&&m==512)
eigS = eigs(T,W);
Umax = max(abs(eigS))
Umin = 0.00006011
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

% disp('cholesky分解：')
% [IT,res]=cholesky(W,T,n,alpha)
% 
% 
% disp('R重排的cholesky分解：')
% [IT,res]=rpcholesky(W,T,n,alpha)

disp('LU分解：')
tic;
[IT,res]=LU(W,T,n,alpha,p,q)
toc;

disp('重排的cholesky分解：')
tic;
[IT,res]=cpcholesky(W,T,n,alpha,p,q)
toc;


