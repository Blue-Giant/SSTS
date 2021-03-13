%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric matrices  are 
%  W = (Im) kron (Bm) + (Bm) kron (Im)is positive   T = (1/m*m)*[(Im) kron (Cm) + (Cm) kron (Im)]
%
%      2     -1      0  ......  0    0             m-1   -1     -1  ......   -1    -1
%      -1     2     -1  ......  0    0             -1   m-1     -1  ......   -1    -1
%Bm =  0     -1      2  ......  0    0       Cm =  -1    -1    m-1  ......   -1    -1
%      ...   ...   ...  .      ...  ...            ...   ...   ...  .        ...  ... 
%      ...   ...   ...    .    ...  ...            ...   ...   ...    .      ...  ...
%      0      0    ...         -1    2              0     0    ...           -1    m-1
%its coresspoding complex symmetric linear systems is
%                             AA = W + iT
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SBTSe5�����Ĳ���ѡȡ�б�
%m         10               20               30             40               50                  
%alpha   1.0029         1.0026               
%beta    0.9971         0.9974            
%Umax    0.0041         0.0037         
%Umin   8.5057e-05      2.0950e-05    
%IT        2               2                11              11                10

%alpha(IT)   2.13 (10)   2.07(11)         2.05(11)         2.04(10)           2.04(10)      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�����������Ҫ��һЩ�������÷�
%��һ��6x6�����У���ȡ���ĵ�һ��Ԫ���γ�һ��6ά��������A(i,:)�� ,A(:,j)�з�����
%A(i,:)��ȡA�ĵ�i��
%A(:,j)��ȡA�ĵ�j��
clear;                                      %����ϴ������ɵ��ڴ��еĶ���
clc;                                        %����
m=input('please input  m=');                %��8��m��ֵ=8,16,32,64,128,256
OMEGA=input('please input  OMEGA=');        %����alpha��ֵ
alpha=input('please input  alpha=');        %�������alpha��ֵ���������alpha=0�������������alpha��ֵ
h=1/(m+1);                                      %���񲽳�
n=m*m;                                      %���ɾ���ά��
Bm=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%���ԽǾ���(-1,2-1)                            
W =  kron(speye(m),Bm) + kron(Bm,speye(m));

Cm = m*speye(m) - ones(m);
T = ( 1/m*m )*(  kron(speye(m),Cm) + kron(Cm,speye(m)) );

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
