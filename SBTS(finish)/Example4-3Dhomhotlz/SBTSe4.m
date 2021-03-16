%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex 3-D Helmholtz equations is 
%   -\LapLace U + k*k U +i\sigma2 U = f(x,y,z)  (x,y,z)in \Omega=[0,1]x[0,1]x[0,1]
%    U|F = g(x,y,z),    (x,y,z) in F
%its coresspoding complex symmetric linear systems is
%                             A = W + iT
%  W = (Vm)kron(Im)kron(Im) + (Im)kron(Vm)kron(Im) + (Im)kron(Im)kron(Vm) - k*k*h*h[(Im)kron(Im)kron(Im)]   
%       Vm = tridiag(-1,2-1)
%  T = \sigma2[(Im)kron(Im)kron(Im)] 
%         h=1/m      k= 1, sigma2 =0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SBTS方法的参数选取列表
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
k = 1;
sigma2 = 0.1;
m=input('please input  m=');                %输8入m的值=8,16,32,64,128,256
alpha=input('please input  alpha=');        %输入参数alpha的值，如果输入alpha=0，程序会计算参数alpha的值
h=1/(m+1);                                  %网格步长
n=m*m*m;                                    %生成矩阵维数
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%三对角矩阵(-1,2-1)
Vm=(1/(h*h))*V;                             
K = kron(kron(Vm,speye(m)),speye(m))+  kron(kron(speye(m),Vm),speye(m)) + kron(kron(speye(m),speye(m)),Vm);
IMMM = k*k*h*h*(     kron(kron(speye(m),speye(m)),speye(m))   );

W = K-IMMM;
W = h*h*W;
T = sigma2*(     kron(kron(speye(m),speye(m)),speye(m))   );
T = h*h*T;


if (alpha==0)
S=W\T;
eigS = eig(full(S));
Umax = max(abs(eigS))
Umin = min(abs(eigS))
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

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


