%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      GSOE  method
%The complex symmetric linear systems of equations is  
%                 [-(w^2)M +K +i(wCv + Ch)]x = b
%     W =-(w^2)M +K  T = wCv + Ch    with Cv = 10I(��λ��)     Ch =uK
%     M  = I(��λ��)      K=I\kron Vm + Vm\kron I
%     Vm = (h^-2)*[tridiag(-1,2,-1)]   h=1/m 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GSOR�����Ĳ���ѡȡ�б�
%m         8          16       32        48       64        80       128       256  
%alpha   0.451      0.455    0.457      0.457    0.457     0.457    0.457     0.457   
%IT        32         26       28        26       25        24        23        23
%RES     8.6656     7.525    7.1347    6.1490    5.9843   7.8108    9.9998    8.8207
%CPU     0.00107   0.00137   0.00724   0.01933   0.0425   0.0848    0.3238    2.2089��ʱ����Ҫ�����ܣ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��һ��6x6�����У���ȡ���ĵ�һ��Ԫ���γ�һ��6ά��������A(i,:)�� ,A(:,j)�з�����
%A(i,:)��ȡA�ĵ�i��
%A(:,j)��ȡA�ĵ�j��
clear;                                      %����ϴ������ɵ��ڴ��еĶ���
clc;                                        %����
m=input('please input  m=');                %����m��ֵ=8,16,32,64,128,256
OMEGA=input('please input  OMEGA=');        %����alpha��ֵ,���alpha=0,��������Ų���
alpha=input('please input  alpha=');        %����alpha��ֵ
h=1/(m+1);                                      %���񲽳�
n=m*m;                                      %���ɵľ���ά��
omega = pi;
mu = 0.02;
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%���ԽǾ���(-1,2,-1)
Vm=(1/(h*h))*V;                             
K=kron(speye(m),Vm)+kron(Vm,speye(m));

Cv = 10*speye(n);
Ch = mu*K;
M = speye(n);
W=-omega*omega*M + K;
W=h*h*W;
T=omega*Cv + Ch;
T=h*h*T;

p=(W-T)*ones(n,1);
q=(W+T)*ones(n,1);


%��Ԥ�������
if  (OMEGA==0&&m~=128&&m~=256&&m~=512)
S=W\T;
eigS = eig(full(S));
umax = max(abs(eigS));
umin = min(abs(eigS));
OMEGA_fenzi = 1- umax*umin+sqrt((1+umin*umin)*(1+umax*umax)) ;
OMEGA_fenmu = umax + umin;
OMEGA=OMEGA_fenzi/OMEGA_fenmu
end

%Ԥ������ϵ��������Ҷ�����
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

disp('LU�ֽ⣺')
[IT,res]=LU(WW,TT,n,alpha,pp,qq)

disp('PCPLU�ֽ⣺')
tic;
[ITPCP,resPCP]=PCPLU(WW,TT,n,alpha,pp,qq);
toc;
ITPCP
resPCP

disp('chol�ֽ⣺')
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

