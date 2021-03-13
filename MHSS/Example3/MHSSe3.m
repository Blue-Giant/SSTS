%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is          
%                           (W+iT)*u = b; 
%     W=10(I\kron Vc + Vc\kron I)+9(e1*em'+em*e1')     T=I\kron Vc + Vc\kron I
%     V = tridiag(-1,2,-1)              Vc = V -  e1*em'+em*e1'   
%     e1 单位向量(1,0,0,0,......,,0)    em 单位向量(0,0,0,0,......,,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  A = (W+iT)x = [(alphaI+W)-(alphaI-iT)]x = b  (1)
%  -iA = (T-iW) = [(alphaI+T)-(alphaI+iW)]x = -ib (2)
% 根据(1) 和 （2）构造了迭代格式
% (alphaI+W)x^(k+0.5) = (alphaI-iT)x^(k)+b
% (alphaI+T)x^(k+1) = (alphaI+iW)x^(k+0.5)-ib
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HSSE3方法的参数选取列表
%m         8        16        32       48       64       80      128        256  
%alpha    3.09     1.61      1.01     0.73     0.53     0.44     0.26       0.13
%IT        36       53        76      100      130      156       246       468 
%res     9.5702    9.4687    9.0949  9.6682   9.7030   9.9911    9.9911    9.8256
%cpu     0.0011    0.0044    0.0332  0.1411   0.4136   0.9552    5.6178    76.8932                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%从一个6x6矩阵中，提取它的第一行元素形成一个6维行向量。A(i,:)行 ,A(:,j)列方法。
%A(i,:)提取A的第i行
%A(:,j)提取A的第j列
clear;                                      %清除上次所生成的内存中的东西
clc;                                        %清屏
m=input('please input  m=');                %输入m的值=8,16,32,64,128,256
alpha=input('please input  alpha=');        %请输入参数alpha的值
h=1/(m+1);                                      %网格步长
n=m*m;                                      %生成矩阵维数
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%三对角矩阵(1-,2-1)
e=speye(m);                                 %离散稀疏存储的单位矩阵
e1=e(:,1);                                  %单位向量(1,0,0,0,......,,0)
em=e(:,m);                                  %单位向量(0,0,0,0,......,,1)
Ve=e1*em'+em*e1';
Vc=V-Ve;

W=10*(kron(e,Vc)+kron(Vc,e)) +9*kron(Ve,e); %对称正定矩阵
% W=h*h*W;   
T=kron(e,V)+kron(V,e);                      %对称矩阵
% T=h*h*T;

%矩阵分块命令mat2cell
% mat2cell(W,())

% disp('cholesky分解：')
% [IT,res]=cholesky(W,T,n,alpha)
% 
% disp('cpcholesky分解：')
% [IT,res]=cpcholesky(W,T,n,alpha)
% 
% disp('rpcholesky分解：')
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



