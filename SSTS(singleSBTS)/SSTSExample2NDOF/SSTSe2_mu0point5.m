%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is  
%                 [-(w^2)M +K +i(wCv + Ch)]x = b
%     W =-(w^2)M +K  T = wCv + Ch    with Cv = 10I(��λ��)     Ch =uK
%     M  = I(��λ��)      K=I\kron Vm + Vm\kron I
%     Vm = (h^-2)*[tridiag(-1,2,-1)]   h=1/m        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%preconditioned single SBTS�����Ĳ���ѡȡ�б�
%m           8         16        32         64        128        256        512
%alpha     1.1923    1.1897    1.1889     1.1887     1.1887    1.1887      1.1887
%eta_max   0.6197    0.6157    0.6145     0.6142     0.6142    0.6184      0.6142
%eta_min   0.0250    0.0163    0.0138     0.0131     0.0130    0.01295     0.0129
%IT          8          8         7          7          6        6           6
%res       6.0991    2.5282    6.0804     2.2302     5.9391     3.4222     2.9319
%----------------ȡС�������λ��Ч��ֵ--------------------------------------------------------------------
%alpha(IT) 1.19(8)   1.19(8)    1.19(7)   1.19(7)    1.19(6)    1.19(6)    1.19(6)
%res       6.6405    2.4994     5.9010    2.1556     5.9453     3.3356     2.8173
%cpu(s)    0.0001    0.0003     0.0015    0.0108     0.0762     0.5237     4.1122
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��һ��6x6�����У���ȡ���ĵ�һ��Ԫ���γ�һ��6ά��������A(i,:)�� ,A(:,j)�з�����
%A(i,:)��ȡA�ĵ�i��
%A(:,j)��ȡA�ĵ�j��
clear;                                      %����ϴ������ɵ��ڴ��еĶ���
clc;                                        %����
m=input('please input  m=');                %����m��ֵ=8,16,32,64,128,256,512
alpha=input('please input  alpha=');        %���������ֵ
h=1/(m+1);                                  %���񲽳�
n=m*m;                                      %���ɵľ���ά��
omega = pi;
mu = 0.5;
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%���ԽǾ���(-1,2,-1)
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
SS=WW\TT;
eigSS = eig(full(SS));
Etamax = max(abs(eigSS))
Etamin = min(abs(eigSS))
temp = 2+ Etamax*Etamax + Etamin*Etamin;
alpha = (temp)/2
end

if (alpha==0&&m==256)
Etamax = 0.6142
Etamin = 0.01295
temp = 2+ Etamax*Etamax + Etamin*Etamin;
alpha = (temp)/2
end
if (alpha==0&&m==512)
Etamax = 0.6142
Etamin = 0.0129
temp = 2+ Etamax*Etamax + Etamin*Etamin;
alpha = (temp)/2
end

disp('LU�ֽ⣺')
tic;
[IT,res]=LU(WW,TT,n,alpha,pp,qq);
toc;
IT
res

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


