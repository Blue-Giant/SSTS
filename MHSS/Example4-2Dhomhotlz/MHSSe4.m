%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex Helmholtz equations is 
%                   -\LapLace U + \sigma1 U +i\sigma2 U = f
%its coresspoding complex symmetric linear systems is
%                    ((K + \sigam1*I) +i\sigam2*I)*u = b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  A = (W+iT)x = [(alphaI+W)-(alphaI-iT)]x = b  (1)
%  -iA = (T-iW) = [(alphaI+T)-(alphaI+iW)]x = -ib (2)
% 根据(1) 和 （2）构造了迭代格式
% (alphaI+W)x^(k+0.5) = (alphaI-iT)x^(k)+b
% (alphaI+T)x^(k+1) = (alphaI+iW)x^(k+0.5)-ib
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MHSS方法的参数选取列表
%m          8         16        32       64        128       256  
%alpha    2.25       0.37      0.09     0.021     0.005     0.002 
%IT         23        30        36       39        40        41
%RES     6.1193     9.5640    9.2337   9.0305    9.7788    9.3253
%cpu     0.0006     0.0018    0.0079   0.0573    0.4307    2.7208
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%构造矩阵所需要的一些操作及用法
%从一个6x6矩阵中，提取它的第一行元素形成一个6维行向量。A(i,:)行 ,A(:,j)列方法。
%A(i,:)提取A的第i行
%A(:,j)提取A的第j列
clear;                                      %清除上次所生成的内存中的东西
clc;                                        %清屏
sigma1=100;
sigma2=100;
m=input('please input  m=');                %输8入m的值=8,16,32,64,128,256
alpha=input('please input  alpha=');        %输入参数alpha的值
h=1/(m+1);                                      %网格步长
n=m*m;                                      %生成矩阵维数
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%三对角矩阵(-1,2-1)
Vm=(1/(h*h))*V;                             
K=kron(speye(m),Vm)+kron(Vm,speye(m));


W=K+sigma1*speye(n);
W=h*h*W;
T=sigma2*speye(n);
T=h*h*T;

%矩阵分块命令mat2cell
% mat2cell(W,())

% disp('cholesky分解：')
% [IT,res]=cholesky(W,T,n,alpha)
% 
% disp('重排的cholesky分解：')
% [IT,res]=cpcholesky(W,T,n,alpha)
% 
% disp('R重排的cholesky分解：')
% [IT,res]=rpcholesky(W,T,n,alpha)

disp('LU分解：')
[IT,res]=LU(W,T,n,alpha)


A = W+1i*T;
b = p+1i*q;
tol = 1e-6;
disp('matlab system GMRES')
tic;
MP = (alpha*speye(n)+W)*(alpha*speye(n)+T);
gmres(A,b,20,tol,250,MP);  
toc;



