%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is          
%                           (W+iT)*u = b; 
%     W=10(I\kron Vc + Vc\kron I)+9(e1*em'+em*e1')     T=I\kron Vc + Vc\kron I
%     V = tridiag(-1,2,-1)   Vc = V -  e1*em'+em*e1'   
%     e1 单位向量(1,0,0,0,......,,0)    em 单位向量(0,0,0,0,......,,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PSBTS方法的参数选取列表
%m         8          16        32          64         128       256  
%OMEGA   4.5036     3.002     1.9783      1.4366      1.1813
%alpha   1.1278     1.2339    1.4284      1.6760      1.9047      
%beta    0.8982     0.8407    0.7693      0.7126      0.6780     
%IT        3          4         5            7           8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%构造矩阵所需要的一些操作及用法
%从一个6x6矩阵中，提取它的第一行元素形成一个6维行向量。A(i,:)行 ,A(:,j)列方法。
%A(i,:)提取A的第i行
%A(:,j)提取A的第j列
clear;                                      %清除上次所生成的内存中的东西
clc;                                        %清屏
m=input('please input  m=');                %输入m的值=8,16,32,64,128,256
OMEGA=input('please input  OMEGA=');        %输入参数OMEGA的值,如果输入OMEGA=0，程序会计算参数OMEGA的值
alpha=input('please input  alpha=');        %输入参数alpha的值，如果输入alpha=0，程序会计算参数alpha的值
h=1/m;                                      %网格步长
n=m*m;                                      %生成矩阵维数
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%三对角矩阵(-1,2-1)

e=speye(m);                                 %离散稀疏存储的单位矩阵
e1=e(:,1);                                  %单位向量(1,0,0,...,0)
em=e(:,m);                                  %单位向量(0,0,0,...,1)
Ve=e1*em'+em*e1';
Vc=V-Ve;


W=10*(kron(e,Vc)+kron(Vc,e)) +9*kron(Ve,e); %对称正定矩阵W
W=h*h*W;
T=kron(e,V)+kron(V,e);
T=h*h*T;

%求预处理参数
if (OMEGA==0)
S=W\T;
eigS = eig(full(S));                       %求所有的特征值
umax = max(abs(eigS));                     %特征值绝对值最大数值
umin = min(abs(eigS));                     %特征值绝对值最小数值
OMEGA_fenzi = 1- umax*umin+sqrt((1+umin*umin)*(1+umax*umax)) ;
OMEGA_fenmu = umax + umin;
disp('OMEGA=')                             %在屏幕上显示内容
OMEGA=OMEGA_fenzi/OMEGA_fenmu
end

%预处理后的系数矩阵和右端向量
WW=OMEGA*W+T;
TT=OMEGA*T-W;

%求解正定矩阵W逆的运算，如下
% [Wl,Wu]=lu(W);
% full(Wl);
% full(Wu);
% E=speye(n);
% X=zeros(n);
% Y=zeros(n);
% INV=zeros(n);
% tic;
% for j = 1:n
%     I=E(:,j);
%     full(I);
%      Y=Wl\I;
%      X=Wu\Y;
%     INV(:,j)=X;    
% end
% toc;
% S=INV*T;

%矩阵分块命令mat2cell
% mat2cell(W,())

%最优的迭代参数
if (alpha==0)
SS=WW\TT;
eigSS = eig(full(SS));
Umax = max(abs(eigSS));
Umin = min(abs(eigSS));
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end


disp('LU分解：')
[IT,res]=LU(WW,TT,n,alpha)


