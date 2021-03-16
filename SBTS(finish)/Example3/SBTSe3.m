%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is          
%                           (W+iT)*u = b; 
%     W=10(I\kron Vc + Vc\kron I)+9(e1*em'+em*e1')     T=I\kron Vc + Vc\kron I
%     V = tridiag(-1,2,-1)   Vc = V -  e1*em'+em*e1'   
%     e1 单位向量(1,0,0,0,......,,0)    em 单位向量(0,0,0,0,......,,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SBTS方法的参数选取列表
%m            8          16         32         64         128        256  
%alpha      1.3747     1.7470     2.8821     6.8786     22.1690     82.3060
%beta       0.7858     0.7005     0.6050     0.5392     0.5115      0.5031   
%Umax       0.3962     0.6667     1.2183     2.3270     4.5473      8.9892
%Umin       0.0599     0.0551     0.0526     0.0513     0.0507      0.0500
%IT            5         8          16         43        146         560
%RES        9.5624     5.6641     6.5896     8.4509     9.5010      9.7746
%-----------------------------取两位有效数字-------------------------------
%alpha(IT)  1.37(5)    1.75(8)    2.88(16)   6.88(43)   22.17(146)     48
%res        9.4695     5.8455     6.5230     8.4562     9.4762
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%从一个6x6矩阵中，提取它的第一行元素形成一个6维行向量。A(i,:)行 ,A(:,j)列方法。
%A(i,:)提取A的第i行
%A(:,j)提取A的第j列
clear;                                      %清除上次所生成的内存中的东西
clc;                                        %清屏
m=input('please input  m=');                %输8入m的值=8,16,32,64,128,256
alpha=input('please input  alpha=');        %请输入参数alpha的值,如果alpha=0,程序会自动计算alpha的值
h=1/m;                                      %网格步长
n=m*m;                                      %生成矩阵维数
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%三对角矩阵(1-,2-1)
e=speye(m);                                 %离散稀疏存储的单位矩阵
e1=e(:,1);                                  %单位向量(1,0,0,0,......,,0)
em=e(:,m);                                  %单位向量(0,0,0,0,......,,1)
Ve=e1*em'+em*e1';
Vc=V-Ve;

W=10*(kron(e,Vc)+kron(Vc,e)) +9*kron(Ve,e); %对称正定矩阵
W=h*h*W;                             
T=kron(e,V)+kron(V,e);                      %对称矩阵
T=h*h*T;

%矩阵分块命令mat2cell
% mat2cell(W,())
%求参数值
if (alpha==0&&m~=256)
S=W\T;
eigS = eig(full(S));
Umax = max(abs(eigS))
Umin = min(abs(eigS))
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if (alpha==0&&m==256)
eigS = eigs(T,W);
Umax = max(abs(eigS))
Umin = 0.05
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


