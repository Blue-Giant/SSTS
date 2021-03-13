%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is  
%                 [-(w^2)M +K +i(wCv + Ch)]x = b
%     W =-(w^2)M +K  T = wCv + Ch    with Cv = 10I(单位阵)     Ch =uK
%     M  = I(单位阵)      K=I\kron Vm + Vm\kron I
%     Vm = (h^-2)*[tridiag(-1,2,-1)]   h=1/m        
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Single SBTS方法的参数选取列表
%m             8        16         32          64         128        256        512
%alpha      1.3583    1.4337     1.4562      1.4629      1.4651     1.4672      1.4680
%eta_max    0.8434    0.9312     0.9534      0.9589      0.9603     0.9618      0.9632
%eta_min    0.0729    0.0174     0.0592      0.0792      0.0889     0.0899      0.0912
%IT           9         11         12          13         13          13          13
%res        4.3033    4.8386     7.4969      4.3313      5.4410     5.5531      1.4558
%-----------------------------------------取小数点后两位有效数值-----------------------------------------------
%alpha(IT)  1.36(9)   1.43(11)   1.46(12)    1.46(13)    1.47(13)   1.47(13)    1.47(13)
%res        4.2322    5.4747     6.5982      4.7990      4.5761     5.0436       5.1911
%cpu(s)     0.0001    0.0007     0.0025     0.0187      0.1619     1.1785       8.6784 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%从一个6x6矩阵中，提取它的第一行元素形成一个6维行向量。A(i,:)行 ,A(:,j)列方法。
%A(i,:)提取A的第i行
%A(:,j)提取A的第j列
clear;                                      %清除上次所生成的内存中的东西
clc;                                        %清屏
m=input('please input  m=');                %输入m的值=8,16,32,64,128,256,512
OMEGA=input('please input  OMEGA=');        %输入参数OMEGA的值,如果输入OMEGA=0，程序会计算参数OMEGA的值
alpha=input('please input  alpha=');        %输入参数的值
h=1/(m+1);                                      %网格步长
n=m*m;                                      %生成的矩阵维数
omega = pi;
mu = 0.02;
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

%求预处理参数
if (OMEGA==0)
S=W\T;
eigS = eig(full(S));                       %求所有的特征值
umax = max(abs(eigS));                     %特征值绝对值最大数值
umin = min(abs(eigS));                     %特征值绝对值最小数值
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

disp('chol分解：')
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
