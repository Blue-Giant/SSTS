%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is  
%                 [-(w^2)M +K +i(wCv + Ch)]x = b
%     W =-(w^2)M +K  T = wCv + Ch    with Cv = 10I(单位阵)     Ch =uK
%     M  = I(单位阵)      K=I\kron Vm + Vm\kron I
%     Vm = (h^-2)*[tridiag(-1,2,-1)]   h=1/m        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%preconditioned single SBTS方法的参数选取列表
%m            8        16        32         64          128         256         512
%alpha      1.2324   1.2297    1.2291     1.2289       1.2289      1.2228      1.2228
%eta_max    0.6810   0.6778    0.6769     0.6766       0.6766      0.6766      0.6766
%eta_min    0.0323   0.0090    0.0024    6.1089e-04   1.5507e-04 0.3934e-04   9.8400e-06
%IT           9        9         9          9            9           8            8
%res        3.6367   2.2842    3.4827     1.9146       1.9070     2.0514       8.6009
%----------------取小数点后两位有效数值--------------------------------------------------------------------
%alpha(IT) 1.23(9)   1.23(9)    1.23(9)   1.23(9)     1.23(9)      1.22(8)      1.22(8)
%res       3.9385    2.2858     2.0135    1.9808      1.9770       7.9086       7.9065
%cpu(s)    0.0001    0.0004     0.0018    0.0123      0.1071       0.7558       5.0664
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%从一个6x6矩阵中，提取它的第一行元素形成一个6维行向量。A(i,:)行 ,A(:,j)列方法。
%A(i,:)提取A的第i行
%A(:,j)提取A的第j列
clear;                                      %清除上次所生成的内存中的东西
clc;                                        %清屏
m=input('please input  m=');                %输入m的值=8,16,32,64,128,256,512
alpha=input('please input  alpha=');        %输入参数的值
h=1/(m+1);                                      %网格步长
n=m*m;                                      %生成的矩阵维数
omega = pi;
mu = 1;
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

WW = W+T;
TT = T-W;
pp = p+q;
qq = q-p;

if (alpha==0&&m~=128&&m~=256&&m~=512)
SS=WW\TT;
eigSS = eig(full(SS));
Etamax = max(abs(eigSS))
Etamin = min(abs(eigSS))
temp = 2+ Etamax*Etamax + Etamin*Etamin;
alpha = (temp)/2
end

if (alpha==0&&m==128)
% SS=WW\TT;
% eigSS = eig(full(SS));
Etamax = 0.6766
Etamin = 1.5507e-04
temp = 2+ Etamax*Etamax + Etamin*Etamin;
alpha = (temp)/2
end

if (alpha==0&&m==256)
Etamax = 0.6676
Etamin = 0.3934e-04
temp = 2+ Etamax*Etamax + Etamin*Etamin;
alpha = (temp)/2
end
if (alpha==0&&m==512)
Etamax = 0.6676
Etamin = 0.0984e-04
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

disp('chol分解：')
tic;
[ITchol,reschol]=cpcholesky(WW,TT,n,alpha,pp,qq);
toc;
ITchol
reschol


