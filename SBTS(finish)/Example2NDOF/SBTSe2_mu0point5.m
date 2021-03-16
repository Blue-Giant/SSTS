%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is  
%                 [-(w^2)M +K +i(wCv + Ch)]x = b
%     W =-(w^2)M +K  T = wCv + Ch    with Cv = 10I(单位阵)     Ch =uK
%     M  = I(单位阵)      K=I\kron Vm + Vm\kron I
%     Vm = (h^-2)*[tridiag(-1,2,-1)]   h=1/(m+1)        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SBTS方法的参数选取列表
%m           8          16         32         64          128         256         512
%alpha     19.9395    19.4279    19.2858    19.2481      19.2383    19.2358     19.2352
%beta      0.5129     0.5132     0.5133     0.5133       0.5133     0.5133      0.5133
%Umax      4.2591     4.2042     4.1887     4.1845       4.1835     4.1832      4.1831 
%Umin      0.5588     0.5159     0.5042     0.5011       0.5003     0.5001      0.5000
%IT          104       100        100        100           101        101        101
%res       9.0280     9.7181     9.0472     9.8685       8.9729     9.0915      9.1250
%---------------------------------------------------------------------------------------------------------------------
%对比MHSS
% IT                    51         85        150            279       529
%----------------取小数点后两位有效数值--------------------------------------------------------------------
%alpha(IT) 19.94(104) 19.43(100) 19.29(100)  19.25(100)  19.24(101)  19.24(103)  19.24(101)
%res       8.9833      9.5806     9.0242     9.8796       8.9840      9.1195     9.1578
%CPU(s)    0.0026      0.0061     0.0318     0.2267       1.9161      13.7038    102.2811
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
mu = 0.5;
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
Umax = max(abs(eigS))
Umin = min(abs(eigS))
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if(alpha==0&&m==128)
S=W\T;
eigS = eig(full(S));
Umax = max(abs(eigS))
Umin = min(abs(eigS))
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if(alpha==0&&m==256)
eigS = eigs(T,W);
Umax = max(abs(eigS))
Umin = 0.5001
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if(alpha==0&&m==512)
eigS = eigs(T,W);
Umax = max(abs(eigS))
Umin = 0.5
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

disp('LU分解：')
[IT,res]=LU(W,T,n,alpha,p,q)


