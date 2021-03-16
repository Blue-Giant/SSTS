%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex Helmholtz equations is 
%                   -\LapLace U + \sigma1 U +i\sigma2 U = f
%its coresspoding complex symmetric linear systems is
%                   ((K + \sigam1*I) +i\sigam2*I)*u = b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PSBTS�����Ĳ���ѡȡ�б�
%m         8        16         32          64         128       256 
%OMEGA   2.1243    2.5320    2.6859      2.7338      2.7488     
%alpha   1.2420    1.3069    1.3259      1.3301      1.3308       
%beta    0.8369    0.8098    0.8027      0.8012      0.8009    
%IT         5        5         5            5           5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�����������Ҫ��һЩ�������÷�
%��һ��6x6�����У���ȡ���ĵ�һ��Ԫ���γ�һ��6ά��������A(i,:)�� ,A(:,j)�з�����
%A(i,:)��ȡA�ĵ�i��
%A(:,j)��ȡA�ĵ�j��
clear;                                      %����ϴ������ɵ��ڴ��еĶ���
clc;                                        %����
sigma1=100; 
sigma2=100;
m=input('please input  m= ');               %��8��m��ֵ=8,16,32,64,128,256
OMEGA=input('please input  OMEGA= ');       %�������alpha��ֵ���������OMEGA=0�������������OMEGA��ֵ
alpha=input('please input  alpha= ');       %�������alpha��ֵ���������alpha=0�������������alpha��ֵ
h=1/m;                                      %���񲽳�
n=m*m;                                      %���ɾ���ά��
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%���ԽǾ���(-1,2-1)
Vm=(1/(h*h))*V;                             
K=kron(speye(m),Vm)+kron(Vm,speye(m));


W=K+sigma1*speye(n);
W=h*h*W;
T=sigma2*speye(n);
T=h*h*T;

%��Ԥ�������
if (OMEGA==0)
S=W\T;
eigS = eig(full(S));                       %�����е�����ֵ
umax = max(abs(eigS));                     %����ֵ����ֵ�����ֵ
umin = min(abs(eigS));                     %����ֵ����ֵ��С��ֵ
OMEGA_fenzi = 1- umax*umin+sqrt((1+umin*umin)*(1+umax*umax)) ;
OMEGA_fenmu = umax + umin;
disp('OMEGA=')                             %����Ļ����ʾ����
OMEGA=OMEGA_fenzi/OMEGA_fenmu
end

%Ԥ������ϵ��������Ҷ�����
WW=OMEGA*W+T;
TT=OMEGA*T-W;

%����ֿ�����mat2cell
% mat2cell(W,())

%���ŵĵ�������
if (alpha==0)
SS=WW\TT;
eigSS = eig(full(SS));
Umax = max(abs(eigSS));
Umin = min(abs(eigSS));
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

disp('cholesky�ֽ⣺')
[IT,res]=cholesky(WW,TT,n,alpha)

disp('���ŵ�cholesky�ֽ⣺')
[IT,res]=cpcholesky(WW,TT,n,alpha)

disp('R���ŵ�cholesky�ֽ⣺')
[IT,res]=rpcholesky(WW,TT,n,alpha)

disp('LU�ֽ⣺')
[IT,res]=LU(WW,TT,n,alpha)

