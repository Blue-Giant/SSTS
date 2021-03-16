%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is  
%                 (W+iT)x = b
%     W =K+[(3-sqrt(3))/tau]I         T =K+[(3+sqrt(3))/tau]I
%          I(单位阵)                 K=I\kron Vm + Vm\kron I
%     Vm = (h^-2)*[tridiag(-1,2,-1)]   h=1/(m+1)        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SBTS方法的参数选取列表
%m             8           16          32           64          128          256         512
%alpha       6.5880      8.4153      10.6629      12.7600      14.3078      15.2745    15.8198
%beta        0.5411      0.5316      0.5246       0.5204       0.5181       0.5169     0.5163
%Umax        2.0073      2.4280      2.8568       3.2042       3.4379       3.5760     3.6517
%Umin        1.0487      1.0255      1.0131       1.0066       1.0035       1.0018     1.0010
%IT            16          23          31            39          45           48         51
%RES         7.4537      5.8411      6.8003       7.2331       7.6323       9.9263     7.8649
%----------------取小数点后两位有效数值----------------------------------------------
%alpha(IT)  6.59(16)     8.42(23)    10.66(31)    12.76(39)   14.31(45)     15.27(48)  15.82(51)
%res        7.3430       5.7071      6.8119       7.2328      7.6502        9.8802     7.8669
%Cpu(s)     0.0017       0.0045      0.0247       0.2185      1.9279        13.989     108.6717
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%从一个6x6矩阵中，提取它的第一行元素形成一个6维行向量。A(i,:)行 ,A(:,j)列方法。
%A(i,:)提取A的第i行
%A(:,j)提取A的第j列
clear;                                      %清除上次所生成的内存中的东西
clc;                                        %清屏
m=input('please input  m=');                %输入m的值=8,16,32,64,128,256
alpha=input('please input  alpha=');        %输入参数的值
h=1/(m+1);                                  %网格步长
n=m*m;                                      %生成的矩阵维数
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%三对角矩阵(-1,2,-1)
Vm=(1/(h*h))*V;                             
K=kron(speye(m),Vm)+kron(Vm,speye(m));

W=K+((3-sqrt(3))/h)*speye(n);
W=h*h*W;
T=K+((3+sqrt(3))/h)*speye(n);
T=h*h*T;

c=sparse(n,1);
for k=1:n
    c(k)=k/(h*(k+1)*(k+1));   
end

p=h*h*c;
q=-h*h*c;

if (alpha == 0&&m~=128&&m~=256&&m~=512)
disp('sassjd');
S = W\T;
eigS = eig(full(S));
Umax = max(abs(eigS));
Umin = min(abs(eigS));
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if (alpha==0&&m==128)
Umax = max(eigs(T,W));
Umin = 1.0035;
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if (alpha==0&&m==256)
Umax = max(eigs(T,W));
Umin = 1.0018;
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end


if (alpha==0&&m==512)
eigS = eigs(T,W);
Umax = max(abs(eigS));
Umin = 1.0010;
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


disp('******************* mygmres  ***************')
AA = [W -T;T  W];
bb = [p;q];
maxIt = 250;
tol = 1e-6;

% M = [W sparse(n,n); T alpha*W]*[W\speye(n) sparse(n,n); sparse(n,n) (1/(2*alpha-1))*(W\speye(n))]*[W -T;sparse(n,n) alpha*W];
% tic;
% mygmres(AA,bb,10,tol,250,M);
% toc;

disp('******************* SSTSgmres  ***************')
tic;
SBTSgmres(AA,bb,10,tol,250,W,T,alpha);
toc;


