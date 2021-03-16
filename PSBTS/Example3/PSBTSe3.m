%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is          
%                           (W+iT)*u = b; 
%     W=10(I\kron Vc + Vc\kron I)+9(e1*em'+em*e1')     T=I\kron Vc + Vc\kron I
%     V = tridiag(-1,2,-1)   Vc = V -  e1*em'+em*e1'   
%     e1 ��λ����(1,0,0,0,......,,0)    em ��λ����(0,0,0,0,......,,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PSBTS�����Ĳ���ѡȡ�б�
%m         8          16        32          64         128       256  
%OMEGA   4.5036     3.002     1.9783      1.4366      1.1813
%alpha   1.1278     1.2339    1.4284      1.6760      1.9047      
%beta    0.8982     0.8407    0.7693      0.7126      0.6780     
%IT        3          4         5            7           8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�����������Ҫ��һЩ�������÷�
%��һ��6x6�����У���ȡ���ĵ�һ��Ԫ���γ�һ��6ά��������A(i,:)�� ,A(:,j)�з�����
%A(i,:)��ȡA�ĵ�i��
%A(:,j)��ȡA�ĵ�j��
clear;                                      %����ϴ������ɵ��ڴ��еĶ���
clc;                                        %����
m=input('please input  m=');                %����m��ֵ=8,16,32,64,128,256
OMEGA=input('please input  OMEGA=');        %�������OMEGA��ֵ,�������OMEGA=0�������������OMEGA��ֵ
alpha=input('please input  alpha=');        %�������alpha��ֵ���������alpha=0�������������alpha��ֵ
h=1/m;                                      %���񲽳�
n=m*m;                                      %���ɾ���ά��
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%���ԽǾ���(-1,2-1)

e=speye(m);                                 %��ɢϡ��洢�ĵ�λ����
e1=e(:,1);                                  %��λ����(1,0,0,...,0)
em=e(:,m);                                  %��λ����(0,0,0,...,1)
Ve=e1*em'+em*e1';
Vc=V-Ve;


W=10*(kron(e,Vc)+kron(Vc,e)) +9*kron(Ve,e); %�Գ���������W
W=h*h*W;
T=kron(e,V)+kron(V,e);
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

%�����������W������㣬����
% [Wl,Wu]=lu(W);
% full(Wl);
% full(Wu);
% E=speye(n);
% X=zeros(n);
% Y=zeros(n);
% INV=zeros(n);
% tic;
% for j = 1:n
%     I=E(:,j);
%     full(I);
%      Y=Wl\I;
%      X=Wu\Y;
%     INV(:,j)=X;    
% end
% toc;
% S=INV*T;

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


disp('LU�ֽ⣺')
[IT,res]=LU(WW,TT,n,alpha)


