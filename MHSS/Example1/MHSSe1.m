%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is  
%                 (W+iT)x = b
%     W =K+[(3-sqrt(3))/tau]I         T =K+[(3+sqrt(3))/tau]I
%          I(��λ��)                 K=I\kron Vm + Vm\kron I
%     Vm = (h^-2)*[tridiag(-1,2,-1)]   h=1/m        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  A = (W+iT)x = [(alphaI+W)-(alphaI-iT)]x = b  (1)
%  -iA = (T-iW) = [(alphaI+T)-(alphaI+iW)]x = -ib (2)
% ����(1) �� ��2�������˵�����ʽ
% (alphaI+W)x^(k+0.5) = (alphaI-iT)x^(k)+b
% (alphaI+T)x^(k+1) = (alphaI+iW)x^(k+0.5)-ib
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MHSS�����Ĳ���ѡȡ�б�
%m          8         16       32      48        64      80     128      256      512
%alpha     1.66      1.06     0.75    0.66      0.54    0.51   0.40      0.30     0.21
%IT         30        40       54      64        73      80     98       133       181
%RES      6.5749    9.6723   9.6109   8.5006   9.4109  8.8040  9.3469   9.9918    9.7394
%Cpu(s)   0.0016    0.0065   0.0444   0.1791   0.4998  1.0952  5.4905   54.5122  544.3006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%��һ��6x6�����У���ȡ���ĵ�һ��Ԫ���γ�һ��6ά��������A(i,:)�� ,A(:,j)�з�����
%A(i,:)��ȡA�ĵ�i��
%A(:,j)��ȡA�ĵ�j��
clear;                                      %����ϴ������ɵ��ڴ��еĶ���
clc;                                        %����
m=input('please input  m=');                %����m��ֵ=8,16,32,64,128,256
alpha=input('please input  alpha=');        %�������alpha��ֵ
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

b = p+1i*q;


disp('LU�ֽ⣺')
tic;
[IT,res]=LU(W,T,n,alpha,b);
toc;
IT
res

disp('PCPLU�ֽ⣺')
tic;
[ITCPC,resCPC]=PCPLU(W,T,n,alpha,b);
toc;
ITCPC
resCPC

disp('���ŵ�Pcholesky�ֽ⣺')
tic;
[ITchol,reschol]=cpcholesky(W,T,n,alpha,b);
toc;
ITchol
reschol

A = W+1i*T;
b = p+1i*q;
tol = 1e-6;
disp('matlab system GMRES')
tic;
MP = (alpha*speye(n)+W)*(alpha*speye(n)+T);
gmres(A,b,20,tol,250,MP);  
toc;