%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex Helmholtz equations is 
%                   -\LapLace U + \sigma1 U +i\sigma2 U = f
%its coresspoding complex symmetric linear systems is
%                   ((K + \sigam1*I) +i\sigam2*I)*u = b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%  %%%%      %%%%     %%%%   %%%%    %%%%   %%%%   %%%%  %%%%[\
%  % WW    -TT  %     % WW    O %    %  O    O %   % O   TT %
%  %            %  =  %         % -- %         % - %        % = D-L-U
%  % TT     WW  %     % O    WW %    % -TT   O %   % O    O %
%  %%%%      %%%%     %%%%   %%%%    %%%%   %%%%   %%%%  %%%%
% M = (1/alpha)*(D-alpha*L) and is used a preconditioner for (1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PSBTS方法的参数选取列表
%m         8        16         32          64         128       256 
%OMEGA   2.1243    2.5320    2.6859      2.7338      2.7488     
%alpha   1.2420    1.3069    1.3259      1.3301      1.3308       
%beta    0.8369    0.8098    0.8027      0.8012      0.8009    
%IT         5        5         5            5           5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%构造矩阵所需要的一些操作及用法
%从一个6x6矩阵中，提取它的第一行元素形成一个6维行向量。A(i,:)行 ,A(:,j)列方法。
%A(i,:)提取A的第i行
%A(:,j)提取A的第j列
clear;                                      %清除上次所生成的内存中的东西
clc;                                        %清屏
sigma1=100; 
sigma2=100;
m=input('please input  m= ');               %输8入m的值=8,16,32,64,128,256
OMEGA=input('please input  OMEGA= ');       %输入参数alpha的值，如果输入OMEGA=0，程序会计算参数OMEGA的值
alpha=input('please input  alpha= ');       %输入参数alpha的值，如果输入alpha=0，程序会计算参数alpha的值
h=1/m;                                      %网格步长
n=m*m;                                      %生成矩阵维数
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%三对角矩阵(-1,2-1)
Vm=(1/(h*h))*V;                             
K=kron(speye(m),Vm)+kron(Vm,speye(m));

W=K+sigma1*speye(n);
W=h*h*W;
T=sigma2*speye(n);
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
SS = WW\TT;
rhoSS = max( eigs(full(SS)) );
temp = 1+sqrt(1+rhoSS*rhoSS);
alpha = 2/temp
end

if (alpha==0&&m==128)
rhoSS = max( eigs(TT,WW) );
temp = 1+sqrt(1+rhoSS*rhoSS);
alpha = 2/temp
end

disp('LU分解：')
[IT,res]=LU(WW,TT,n,alpha,pp,qq)

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

disp('******************* gmres  ***************')
AA = [WW -TT;TT  WW];
% AA = [W -T;T  W];
M = [WW sparse(n,n); alpha*TT WW];
bb = [pp;qq];
maxIt = 250;
tol = 1e-6;
tic;
gmres(AA,bb,10,tol,250,M);
toc;


disp('******************* mygmres  ***************')
tic;
mygmres(AA,bb,10,tol,250,M);
toc;

disp('******************* SSTSgmres  ***************')
tic;
PGSORgmres(AA,bb,10,tol,250,W,T,alpha,OMEGA);
toc;