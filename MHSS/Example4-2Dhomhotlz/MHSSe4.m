%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex Helmholtz equations is 
%                   -\LapLace U + \sigma1 U +i\sigma2 U = f
%its coresspoding complex symmetric linear systems is
%                    ((K + \sigam1*I) +i\sigam2*I)*u = b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  A = (W+iT)x = [(alphaI+W)-(alphaI-iT)]x = b  (1)
%  -iA = (T-iW) = [(alphaI+T)-(alphaI+iW)]x = -ib (2)
% ����(1) �� ��2�������˵�����ʽ
% (alphaI+W)x^(k+0.5) = (alphaI-iT)x^(k)+b
% (alphaI+T)x^(k+1) = (alphaI+iW)x^(k+0.5)-ib
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MHSS�����Ĳ���ѡȡ�б�
%m          8         16        32       64        128       256  
%alpha    2.25       0.37      0.09     0.021     0.005     0.002 
%IT         23        30        36       39        40        41
%RES     6.1193     9.5640    9.2337   9.0305    9.7788    9.3253
%cpu     0.0006     0.0018    0.0079   0.0573    0.4307    2.7208
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�����������Ҫ��һЩ�������÷�
%��һ��6x6�����У���ȡ���ĵ�һ��Ԫ���γ�һ��6ά��������A(i,:)�� ,A(:,j)�з�����
%A(i,:)��ȡA�ĵ�i��
%A(:,j)��ȡA�ĵ�j��
clear;                                      %����ϴ������ɵ��ڴ��еĶ���
clc;                                        %����
sigma1=100;
sigma2=100;
m=input('please input  m=');                %��8��m��ֵ=8,16,32,64,128,256
alpha=input('please input  alpha=');        %�������alpha��ֵ
h=1/(m+1);                                      %���񲽳�
n=m*m;                                      %���ɾ���ά��
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%���ԽǾ���(-1,2-1)
Vm=(1/(h*h))*V;                             
K=kron(speye(m),Vm)+kron(Vm,speye(m));


W=K+sigma1*speye(n);
W=h*h*W;
T=sigma2*speye(n);
T=h*h*T;

%����ֿ�����mat2cell
% mat2cell(W,())

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


A = W+1i*T;
b = p+1i*q;
tol = 1e-6;
disp('matlab system GMRES')
tic;
MP = (alpha*speye(n)+W)*(alpha*speye(n)+T);
gmres(A,b,20,tol,250,MP);  
toc;



