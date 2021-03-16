%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric matrices  are 
%  W = (Im) kron (Bm) + (Bm) kron (Im)is positive   T = (1/m*m)*[(Im) kron (Cm) + (Cm) kron (Im)]
%
%      2     -1      0  ......  0    0             m-1   -1     -1  ......   -1    -1
%      -1     2     -1  ......  0    0             -1   m-1     -1  ......   -1    -1
%Bm =  0     -1      2  ......  0    0       Cm =  -1    -1    m-1  ......   -1    -1
%      ...   ...   ...  .      ...  ...            ...   ...   ...  .        ...  ... 
%      ...   ...   ...    .    ...  ...            ...   ...   ...    .      ...  ...
%      0      0    ...         -1    2              0     0    ...           -1    m-1
%its coresspoding complex symmetric linear systems is
%                             A = W + iT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SBTSe5方法的参数选取列表
%m         10               20               30             40               50                  
%alpha   1.0029         1.0026               
%beta    0.9971         0.9974            
%Umax    0.0041         0.0037         
%Umin   8.5057e-05      2.0950e-05    
%IT        2               2                11              11                10

%alpha(IT)   2.13 (10)   2.07(11)         2.05(11)         2.04(10)           2.04(10)      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%构造矩阵所需要的一些操作及用法
%从一个6x6矩阵中，提取它的第一行元素形成一个6维行向量。A(i,:)行 ,A(:,j)列方法。
%A(i,:)提取A的第i行
%A(:,j)提取A的第j列
clear;                                      %清除上次所生成的内存中的东西
clc;                                        %清屏
m=input('please input  m=');                %输8入m的值=8,16,32,64,128,256
OMEGA=input('please input  OMEGA=');        %输入alpha的值
alpha=input('please input  alpha=');        %输入参数alpha的值，如果输入alpha=0，程序会计算参数alpha的值
h=1/(m+1);                                      %网格步长
n=m*m;                                      %生成矩阵维数
Bm=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%三对角矩阵(-1,2-1)                            
W =  kron(speye(m),Bm) + kron(Bm,speye(m));

Cm = m*speye(m) - ones(m);
T = ( 1/m*m )*(  kron(speye(m),Cm) + kron(Cm,speye(m)) );

p=(W-T)*ones(n,1);
q=(W+T)*ones(n,1);

%求预处理参数
if  OMEGA==0
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

%最优的迭代参数
if(alpha == 0)
SS=WW\TT;
eigSS = eig(full(SS));
Umax = max(abs(eigSS));
Umin = min(abs(eigSS));
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
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



