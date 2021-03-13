%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is  
%                         (W+iT)x = b
%     W =K+[(3-sqrt(3))/tau]I         T =K+[(3+sqrt(3))/tau]I
%          I(单位阵)                   K=I\kron Vm + Vm\kron I
%     Vm = (h^-2)*[tridiag(-1,2,-1)]   h=1/(m+1）   
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GSOR方法的参数选取列表
%m            8         16         32         64          128        256      512
%Umax       2.0073   2.4280      2.8568     3.2042       3.4379     3.5760   3.6517
%alpha     
%IT             19        23         26        27          26         25        29         
%RES          6.9581    5.5868     5.7229    9.1229      7.2335     5.0466    6.7005
%time         0.00113   0.00333    0.01558   0.05273    0.10777     0.20894   0.9250  (时间重新跑)

%GMRES(10)                2(2)       2(4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%从一个6x6矩阵中，提取它的第一行元素形成一个6维行向量。A(i,:)行 ,A(:,j)列方法。
%A(i,:)提取A的第i行
%A(:,j)提取A的第j列
clear;                                      %清除上次所生成的内存中的东西
clc;                                        %清屏
m=input('please input  m=');                %输入m的值=8,16,32,64,128,256
OMEGA=input('please input  OMEGA=');        %输入alpha的值,如果alpha=0,则计算最优参数
alpha=input('please input  alpha=');        %输入alpha的值,如果alpha=0,则计算最优参数
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