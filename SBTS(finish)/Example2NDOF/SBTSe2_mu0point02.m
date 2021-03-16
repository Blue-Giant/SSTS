%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is  
%                 [-(w^2)M +K +i(wCv + Ch)]x = b
%     W =-(w^2)M +K  T = wCv + Ch    with Cv = 10I(单位阵)     Ch =uK
%     M  = I(单位阵)      K=I\kron Vm + Vm\kron I
%     Vm = (h^-2)*[tridiag(-1,2,-1)]   h=1/(m+1)        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SBTS方法的参数选取列表
%m         8          16         32         64          128         256        512
%alpha   12.3028    11.9861    11.8982    11.8749      11.8688     3.2232   11.8669
%beta    0.5212     0.5218     0.5219     0.5220       0.5220      0.5220    0.5220
%Umax    3.2892     3.2414     3.2279     3.2243       3.2234      3.2232    3.2231
%Umin    0.0711     0.0339     0.0236     0.0209       0.0202      0.0198    0.0196
%IT        82        78          77         77           77          77        77
%res    9.7418      9.3809     9.4342     9.1686       9.1188      9.1079    9.1053
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%从一个6x6矩阵中，提取它的第一行元素形成一个6维行向量。A(i,:)行 ,A(:,j)列方法。
%A(i,:)提取A的第i行
%A(:,j)提取A的第j列
clear;                                      %清除上次所生成的内存中的东西
clc;                                        %清屏
m=input('please input  m=');                %输入m的值=8,16,32,64,128,256
alpha=input('please input  alpha=');        %输入alpha的值
h=1/(m+1);                                      %网格步长
n=m*m;                                      %生成的矩阵维数
omega = pi;
mu = 0.02;
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%三对角矩阵(-1,2,-1)
Vm=(1/(h*h))*V;                             
K=kron(speye(m),Vm)+kron(Vm,speye(m));

M = speye(n);
Cv = 10*speye(n);
Ch = mu*K;
W=-omega*omega*M +K;
W=h*h*W;
T=omega*Cv + Ch;
T=h*h*T;

p=(W-T)*ones(n,1);
q=(W+T)*ones(n,1);

if(alpha==0&&m~=128&m~=256&&m~=512)
S=W\T;
eigS = eig(full(S));
Umax = max(abs(eigS));
Umin = min(abs(eigS));
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if(alpha==0&&m==128)
S=W\T;
eigS = eig(full(S));
Umax = max(abs(eigS));
Umin = min(abs(eigS));
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if(alpha==0&&m==256)
eigS = eigs(T,W);
Umax = max(abs(eigS));
Umin = 0.0198;
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if(alpha==0&&m==512)
eigS = eigs(T,W);
Umax = max(abs(eigS));
Umin = 0.0196;
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end


disp('LU分解：')
tic;
[IT,res]=LU(W,T,n,alpha,p,q);
toc;
IT
res

disp('PCPLU分解：')
tic;
[ITPCP,resPCP]=PCPLU(W,T,n,alpha,p,q);
toc;
ITPCP
resPCP

disp('重排的Pcholesky分解：')
tic;
[ITchol,reschol]=cpcholesky(W,T,n,alpha,p,q);
toc;
ITchol
reschol


