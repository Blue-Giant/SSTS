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
%                             A = W + iT
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
alpha=input('please input  alpha=');        %�������alpha��ֵ���������alpha=0�������������alpha��ֵ
h=1/m;                                      %���񲽳�
n=m*m;                                      %���ɾ���ά��
Bm=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%���ԽǾ���(-1,2-1)                            
W =  kron(speye(m),Bm) + kron(Bm,speye(m));

Cm = m*speye(m) - ones(m);
T = ( 1/m*m )*(  kron(speye(m),Cm) + kron(Cm,speye(m)) );

p=(W-T)*ones(n,1);
q=(W+T)*ones(n,1);

% EW = eig(W)
% ET = eig(T)
rank(full(W))
rank(full(T))
% 
% %����ֿ�����mat2cell
% % mat2cell(W,())

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
[IT,res]=LU(W,T,n,alpha,p,q)


