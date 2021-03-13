%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is  
%                 [-(w^2)M +K +i(wCv + Ch)]x = b
%     W =-(w^2)M +K  T = wCv + Ch    with Cv = 10I(单位阵)     Ch =uK
%     M  = I(单位阵)      K=I\kron Vm + Vm\kron I
%     Vm = (h^-2)*[tridiag(-1,2,-1)]   h=1/m        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%preconditioned single SBTS方法的参数选取列表
%m            8        16        32         64          128         256         512
%alpha      1.3495   1.3433    1.3416     1.3412       1.3410      1.3409      1.3409
%eta_max    0.7587   0.7565    0.7558     0.7557       0.7556      0.7555      0.7555
%eta_min    0.3512   0.3383    0.3346     0.3337       0.3334      0.33334     0.33332
%IT           8        8         9          9            9            9           9
%res        9.0091   8.7970    1.7381     1.8981       1.9542      1.9673      1.9719
%-----------------------------------------取小数点后两位有效数值-----------------------------------------------
%alpha(IT)  1.35(9)  1.34(9)   1.34(9)    1.34(9)      1.34(9)     1.34(9)     1.34(9)
%res        8.9114   8.4375    9.8892     1.8332       1.8946      1.9127      1.9176
%cpu(s)     0.0001   0.0004    0.0018     0.0129       0.1075      0.7922      6.0754
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
mu = 2;
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
Etamax = 0.7556
Etamin = 0.3334
temp = 2+ Etamax*Etamax + Etamin*Etamin;
alpha = (temp)/2
end

if (alpha==0&&m==256)
Etamax = 0.7555
Etamin = 0.33334
temp = 2+ Etamax*Etamax + Etamin*Etamin;
alpha = (temp)/2
end
if (alpha==0&&m==512)
Etamax = 0.7555
Etamin = 0.33332
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


