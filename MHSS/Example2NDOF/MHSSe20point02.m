%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is  
%                 [-(w^2)M +K +i(wCv + Ch)]x = b
%     W =-(w^2)M +K  T = wCv + Ch    with Cv = 10I(单位阵)     Ch =uK
%     M  = I(单位阵)      K=I\kron Vm + Vm\kron I
%     Vm = (h^-2)*[tridiag(-1,2,-1)]   h=1/m        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  A = (W+iT)x = [(alphaI+W)-(alphaI-iT)]x = b  (1)
%  -iA = (T-iW) = [(alphaI+T)-(alphaI+iW)]x = -ib (2)
% 根据(1) 和 （2）构造了迭代格式
% (alphaI+W)x^(k+0.5) = (alphaI-iT)x^(k)+b
% (alphaI+T)x^(k+1) = (alphaI+iW)x^(k+0.5)-ib
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MHSS方法的参数选取列表
%m       8        16      32     48       64      80      128     256       512
%alpha  0.58     0.21    0.08   0.05     0.04    0.03    0.02     0.01     0.005
%IT      29       34      38     43       50      57      81      139       250
%res    6.4417  7.5926  7.4211  8.5992  8.4860  9.8750  9.1367   9.6687    9.8620                
%cpu(s) 0.0008  0.0025  0.0120  0.0418  0.1136  0.2453  1.4540   16.8262   520.1294
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

b=p+1i*q;

disp('LU分解：')
tic;
[IT,res]=LU(W,T,n,alpha,b);
toc;
IT
res

disp('PCPLU分解：')
tic;
[ITCPC,resCPC]=PCPLU(W,T,n,alpha,b);
toc;
ITCPC
resCPC

disp('重排的Pcholesky分解：')
tic;
[ITchol,reschol]=cpcholesky(W,T,n,alpha,b);
toc;
ITchol
reschol

A = W+1i*T;
b = p+1i*q;
tol = 1e-6;
disp('matlab system GMRES')
tic;
MP = (alpha*speye(n)+W)*(alpha*speye(n)+T);
gmres(A,b,20,tol,250,MP);  
toc;


