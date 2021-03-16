%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is  
%                 (W+iT)x = b
%     W =K+[(3-sqrt(3))/tau]I         T =K+[(3+sqrt(3))/tau]I
%          I(��λ��)                 K=I\kron Vm + Vm\kron I
%     Vm = (h^-2)*[tridiag(-1,2,-1)]   h=1/m        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PSBTS�����Ĳ���ѡȡ�б�
%m         8        16         32          64          128       256  
%OMEGA   0.6888   0.6529     0.6224      0.6021      0.5904     0.583
%alpha   1.1249   1.1592     1.1874      1.2506      1.2166      
%beta    0.9001   0.8792     0.8637      0.8543      0.8489  
%IT         4        4         4           4           4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%�����������Ҫ��һЩ�������÷�
%��һ��6x6�����У���ȡ���ĵ�һ��Ԫ���γ�һ��6ά��������A(i,:)�� ,A(:,j)�з�����
%A(i,:)��ȡA�ĵ�i��
%A(:,j)��ȡA�ĵ�j��
clear;                                      %����ϴ������ɵ��ڴ��еĶ���
clc;                                        %����
m = input('please input  m=');                %����m��ֵ=8,16,32,64,128,256
OMEGA = input('please input  OMEGA=');        %�������ֵ
alpha = input('please input  alpha=');        %�������ֵ
h=1/(m+1);                                      %���񲽳�
n=m*m;                                      %���ɵľ���ά��
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%���ԽǾ���(-1,2,-1)
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

%��Ԥ�������
if  OMEGA==0
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

%���ŵĵ�������
if(alpha == 0)
SS=WW\TT;
eigSS = eig(full(SS));
Umax = max(abs(eigSS));
Umin = min(abs(eigSS));
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
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

disp('���ŵ�Pcholesky�ֽ⣺')
tic;
[ITchol,reschol]=cpcholesky(WW,TT,n,alpha,pp,qq);
toc;
ITchol
reschol

disp('******************* gmres  ***************')
AA = [WW -TT;TT  WW];
% AA = [W -T;T  W];

bb = [pp;qq];
maxIt = 250;
tol = 1e-6;

% % M = [W sparse(n,n); T alpha*W]*[W\speye(n) sparse(n,n); sparse(n,n) (1/(2*alpha-1))*(W\speye(n))]*[W -T;sparse(n,n) alpha*W];
% tic;
% gmres(AA,bb,10,tol,250,M);
% toc;

% 
% disp('******************* mygmres  ***************')
% tic;
% mygmres(AA,bb,10,tol,250,M);
% toc;

disp('******************* SSTSgmres  ***************')
tic;
PSBTSgmres(AA,bb,10,tol,250,W,T,alpha,OMEGA);
toc;


