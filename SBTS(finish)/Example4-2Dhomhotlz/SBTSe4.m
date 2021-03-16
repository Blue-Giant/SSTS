%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex Helmholtz equations is 
%                   -\LapLace U + \sigma1 U +i\sigma2 U = f
%its coresspoding complex symmetric linear systems is
%                   ((K + \sigam1*I) +i\sigam2*I)*u = b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SBTS�����Ĳ���ѡȡ�б�
%m            8          16         32         64          128           256  
%alpha      2.0582     2.0372     2.0349     2.0346       2.0346       2.0346 
%beta       0.6604     0.6626     0.6629     0.6629       0.6629       0.6629
%Umax       0.8365     0.8355     0.8353     0.8352       0.8352       0.8352
%Umin       0.1373     0.0418     0.0114     0.0030       7.5070e-04   2.0110e-04
%IT          10          10          10        10            10          10
%RES        7.0581     7.9587     8.6917     9.2844       9.7859        9.4414
%-------------------ȡ��λ��Ч����---------------------------------------------------------------------------------------
%alpha(IT)  2.13(10)   2.07(11)   2.05(11)   2.04(10)     2.04(10)  
%res        9.8162     2.7354     2.6081     9.7813       9.9150 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�����������Ҫ��һЩ�������÷�
%��һ��6x6�����У���ȡ���ĵ�һ��Ԫ���γ�һ��6ά��������A(i,:)�� ,A(:,j)�з�����
%A(i,:)��ȡA�ĵ�i��
%A(:,j)��ȡA�ĵ�j��
clear;                                      %����ϴ������ɵ��ڴ��еĶ���
clc;                                        %����
sigma1=100;
sigma2=100;
m=input('please input  m=');                %��8��m��ֵ=8,16,32,64,128,256
alpha=input('please input  alpha=');        %�������alpha��ֵ���������alpha=0�������������alpha��ֵ
h=1/(m+1);                                      %���񲽳�
n=m*m;                                      %���ɾ���ά��
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%���ԽǾ���(-1,2-1)
Vm=(1/(h*h))*V;                             
K=kron(speye(m),Vm)+kron(Vm,speye(m));


W=K+sigma1*speye(n);
W=h*h*W;
T=sigma2*speye(n);
T=h*h*T;


if(alpha==0&&m~=128&m~=256&&m~=512)
S=W\T;
eigS = eig(full(S));
Umax = max(abs(eigS))
Umin = min(abs(eigS))
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if(alpha==0&&m==128)
S=W\T;
eigS = eig(full(S));
Umax = max(abs(eigS))
Umin = min(abs(eigS))
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if(alpha==0&&m==256)
eigS = eigs(T,W);
Umax = max(abs(eigS))
Umin = 0.0002011
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if(alpha==0&&m==512)
eigS = eigs(T,W);
Umax = max(abs(eigS))
Umin = 0.00006011
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

% disp('cholesky�ֽ⣺')
% [IT,res]=cholesky(W,T,n,alpha)
% 
% disp('���ŵ�cholesky�ֽ⣺')
% [IT,res]=cpcholesky(W,T,n,alpha)
% 
% disp('R���ŵ�cholesky�ֽ⣺')
% [IT,res]=rpcholesky(W,T,n,alpha)

disp('LU�ֽ⣺')
[IT,res]=LU(W,T,n,alpha)


