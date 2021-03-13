%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is          
%                           (W+iT)*u = b; 
%     W=10(I\kron Vc + Vc\kron I)+9(e1*em'+em*e1')     T=I\kron Vc + Vc\kron I
%     V = tridiag(-1,2,-1)   Vc = V -  e1*em'+em*e1'   
%     e1 单位向量(1,0,0,0,......,,0)    em 单位向量(0,0,0,0,......,,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       %%%%        %%%%  %%%%    %%%%     %%%%%%             %%%%%%
%       % omegaI    I  %  % W    -T  %     %omegaW+T   -(omegaT-W) %
%  PA = %              %  %          %  =  %                       % ===
%       % -I    omegaI %  % T     W  %     % omegaT-W    omegaW+T  %
%       %%%%        %%%%  %%%%    %%%%     %%%%%%             %%%%%%

%   %%%%        %%%%  %% %%     %%%%    %%%%  
%   % omegaI    I  %  % p %     % omegap+q %      WW = omegaW+T  TT= omegaT-W
%   %              %  %   %  =  %          % (1)  pp = omegap+q  qq = megaq-p
%   % -I    omegaI %  % q %     % omegaq-p %
%   %%%%        %%%%  %% %%     %%%%    %%%%  

%  %%%%      %%%%     %%%%%%   %%%%%%    %%%%            %%%%
%  % WW    -TT  %     % WW       O  %    %  O            TT %
%  %            %  =  %             % -- %                  % = M-N
%  % TT     WW  %     % TT  alphaWW %    %  O   (alpha-1)WW %
%  %%%%      %%%%     %%%%%%   %%%%%%    %%%%            %%%%
% M is used a preconditioner for (1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SSTS方法的参数选取列表
%m          8        16        32         64        128        256        512
%alpha    1.4869   1.4209    1.4099     1.4080     1.4089    1.4093      1.4095
%eta_max  0.8870   0.8955    0.9000     0.9024     0.9036    0.9042      0.9045
%eta_min  0.4325   0.2000    0.0984     0.0418     0.0369    0.0332      0.0305
%IT         9         10       10         9           9        9           10
%res      3.5983   3.7996    3.6666     5.3198     4.6294    3.7285      6.2978
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%从一个6x6矩阵中，提取它的第一行元素形成一个6维行向量。A(i,:)行 ,A(:,j)列方法。
%A(i,:)提取A的第i行
%A(:,j)提取A的第j列
clear;                                      %清除上次所生成的内存中的东西
clc;                                        %清屏
m=input('please input  m=');                %输8入m的值=8,16,32,64,128,256
OMEGA = input('please input  OMEGA=');
alpha=input('please input  alpha=');        %请输入参数alpha的值,如果alpha=0,程序会自动计算alpha的值
h=1/(m+1);                                      %网格步长
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

p=(W-T)*ones(n,1);
q=(W+T)*ones(n,1);

%求预处理参数
if  (OMEGA==0&&m~=128&&m~=256&&m~=512)
S=W\T;
eigS = eig(full(S));
umax = max(abs(eigS));
umin = min(abs(eigS));
OMEGA_fenzi = 1- umax*umin+sqrt((1+umin*umin)*(1+umax*umax)) ;
OMEGA_fenmu = umax + umin;
OMEGA=OMEGA_fenzi/OMEGA_fenmu
end

%预处理后的系数矩阵和右端向量
WW=OMEGA*W+T;
TT=OMEGA*T-W;

pp=OMEGA*p+q;
qq=OMEGA*q-p;

if (alpha==0&&m~=128&&m~=256&&m~=512)
SS=WW\TT;
eigS = eig(full(SS));
Etamax = max(abs(eigS));
Etamin = min(abs(eigS));
temp = 2+ Etamax*Etamax + Etamin*Etamin;
alpha = (temp)/2
end

if (alpha==0&&m==128)
SS=WW\TT;
eigS = eig(full(SS));
Etamax = max(abs(eigS));
Etamin = min(abs(eigS));
temp = 2+ Etamax*Etamax + Etamin*Etamin;
alpha = (temp)/2
end

if (alpha==0&&m==256)
eigS = eigs(TT,WW);
Etamax = max(abs(eigS));
Etamin = 0.0009;
temp = 2+ Etamax*Etamax + Etamin*Etamin;
alpha = (temp)/2
end

if (alpha==0&&m==512)
eigS = eigs(TT,WW);
Etamax = max(abs(eigS));
Etamin = 0.0005;
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

disp('重排的Pcholesky分解：')
tic;
[ITchol,reschol]=cpcholesky(WW,TT,n,alpha,pp,qq);
toc;
ITchol
reschol


disp('******************* mygmres  ***************')
AA = [WW -TT;TT  WW];
% AA = [W -T;T  W];
M = [WW sparse(n,n); TT alpha*WW];
bb = [pp;qq];
maxIt = 250;
tol = 1e-6;
tic;
mygmres(AA,bb,10,tol,250,M);
toc;

disp('******************* SSTSgmres  ***************')
tic;
SSTSgmres(AA,bb,10,tol,250,W,T,alpha,OMEGA);
toc;




