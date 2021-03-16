%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex 3-D Helmholtz equations is 
%   -\LapLace U + k*k U +i\sigma2 U = f(x,y,z)  (x,y,z)in \Omega=[0,1]x[0,1]x[0,1]
%    U|F = g(x,y,z),    (x,y,z) in F
%its coresspoding complex symmetric linear systems is
%                             A = W + iT
%  W = (Vm)kron(Im)kron(Im) + (Im)kron(Vm)kron(Im) + (Im)kron(Im)kron(Vm) - k*k*h*h[(Im)kron(Im)kron(Im)]   
%       Vm = tridiag(-1,2-1)
%  T = \sigma2[(Im)kron(Im)kron(Im)] 
%         h=1/m      k= 1, sigma2 =0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SBTS�����Ĳ���ѡȡ�б�
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
k = 1;
sigma2 = 0.1;
m=input('please input  m=');                %��8��m��ֵ=8,16,32,64,128,256
alpha=input('please input  alpha=');        %�������alpha��ֵ���������alpha=0�������������alpha��ֵ
h=1/(m+1);                                  %���񲽳�
n=m*m*m;                                    %���ɾ���ά��
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%���ԽǾ���(-1,2-1)
Vm=(1/(h*h))*V;                             
K = kron(kron(Vm,speye(m)),speye(m))+  kron(kron(speye(m),Vm),speye(m)) + kron(kron(speye(m),speye(m)),Vm);
IMMM = k*k*h*h*(     kron(kron(speye(m),speye(m)),speye(m))   );

W = K-IMMM;
W = h*h*W;
T = sigma2*(     kron(kron(speye(m),speye(m)),speye(m))   );
T = h*h*T;


if (alpha==0)
S=W\T;
eigS = eig(full(S));
Umax = max(abs(eigS))
Umin = min(abs(eigS))
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


