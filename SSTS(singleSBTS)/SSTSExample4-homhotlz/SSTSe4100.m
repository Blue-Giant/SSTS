%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex Helmholtz equations is 
%                   -\LapLace U + \sigma1 U +i\sigma2 U = f
%its coresspoding complex symmetric linear systems is
%                   ((K + \sigam1*I) +i\sigam2*I)*u = b
%sigma1=100   sigma2=100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%preconditioned single SBTS方法的参数选取列表
%m          8        16        32         64        128        256        512
%alpha    1.2567   1.4176    1.4800     1.4978     1.5024    
%eta_max  0.7129   0.9103    0.9761     0.9939     0.9985    
%eta_min  0.0717   0.0802    0.0849     0.0873     0.0886     
%IT         9         11       12         13         13        
%res      3.1412   5.4106    3.5795     4.8469     9.0112     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%构造矩阵所需要的一些操作及用法
%从一个6x6矩阵中，提取它的第一行元素形成一个6维行向量。A(i,:)行 ,A(:,j)列方法。
%A(i,:)提取A的第i行
%A(:,j)提取A的第j列
clear;                                      %清除上次所生成的内存中的东西
clc;                                        %清屏
sigma1=-10;
sigma2=10;
m=input('please input  m=');                %输8入m的值=8,16,32,64,128,256
OMEGA = input('please input  OMEGA=');
alpha=input('please input  alpha=');        %输入参数的值
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

%求预处理参数
if  (OMEGA==0&&m~=128&&m~=256&&m~=512)
S=W\T;
eigS = eig(full(S));
umax = max(abs(eigS));
umin = min(abs(eigS));
OMEGA_fenzi = 1- umax*umin+sqrt((1+umin*umin)*(1+umax*umax)) ;
OMEGA_fenmu = umax + umin;
OMEGA=OMEGA_fenzi/OMEGA_fenmu
end

%预处理后的系数矩阵和右端向量
WW=OMEGA*W+T;
TT=OMEGA*T-W;

pp=OMEGA*p+q;
qq=OMEGA*q-p;

if (alpha==0&&m~=128&&m~=256&&m~=512)
SS=WW\TT;
eigS = eig(full(SS));
Etamax = max(abs(eigS));
Etamin = min(abs(eigS));
temp = 2+ Etamax*Etamax + Etamin*Etamin;
alpha = (temp)/2
end

if (alpha==0&&m==128)
SS=WW\TT;
eigS = eig(full(SS));
Etamax = max(abs(eigS));
Etamin = min(abs(eigS));
temp = 2+ Etamax*Etamax + Etamin*Etamin;
alpha = (temp)/2
end

if (alpha==0&&m==256)
eigS = eigs(TT,WW);
Etamax = max(abs(eigS));
Etamin = 0.0009;
temp = 2+ Etamax*Etamax + Etamin*Etamin;
alpha = (temp)/2
end

if (alpha==0&&m==512)
eigS = eigs(TT,WW);
Etamax = max(abs(eigS));
Etamin = 0.0005;
temp = 2+ Etamax*Etamax + Etamin*Etamin;
alpha = (temp)/2
end

disp('LU分解：')
tic;
[IT,res]=LU(WW,TT,n,alpha,pp,qq);
toc;
IT
res

disp('PCPLU分解：')
tic;
[ITPCP,resPCP]=PCPLU(WW,TT,n,alpha,pp,qq);
toc;
ITPCP
resPCP

disp('重排的Pcholesky分解：')
tic;
[ITchol,reschol]=cpcholesky(WW,TT,n,alpha,pp,qq);
toc;
ITchol
reschol


