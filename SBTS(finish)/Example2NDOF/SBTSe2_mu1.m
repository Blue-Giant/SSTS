%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is  
%                 [-(w^2)M +K +i(wCv + Ch)]x = b
%     W =-(w^2)M +K  T = wCv + Ch    with Cv = 10I(��λ��)     Ch =uK
%     M  = I(��λ��)      K=I\kron Vm + Vm\kron I
%     Vm = (h^-2)*[tridiag(-1,2,-1)]   h=1/(m+1)        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SBTS�����Ĳ���ѡȡ�б�
%m           8           16          32           64          128           256          512
%alpha     30.3967     29.6410     29.4310      29.3752      29.3608      29.3571      29.3562
%beta      0.5084      0.5086      0.5086       0.5087       0.5087       0.5087       0.5087
%Umax      5.2695      5.2070      5.1894       5.1847       5.1835       5.1832       5.21831
%Umin      1.0667      1.0181      1.0048       1.0012       1.0003       1.00007      1.00002
%IT          97         95           95           96           96           97           97
%res       8.9626     9.0604       9.2639       9.2565       9.8406       8.6887       8.7372
%----------------------------------------------------------------------------------------------------------------------------
%�Ա�MHSS
% IT                     56          94         168           306          555
%----------------ȡС�������λ��Ч��ֵ--------------------------------------------------------------------
%alpha(IT) 30.40(97)    29.64(95)   29.43(95)   29.38(96)     29.36(96)    29.36(97)   29.36(97)
%res        8.7655      9.1004      9.2676      9.2746        9.8369        8.7014      8.7540
%cpu(s)     0.0024      0.0057      0.0301      0.2133        1.8404        13.1137     101.9332
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��һ��6x6�����У���ȡ���ĵ�һ��Ԫ���γ�һ��6ά��������A(i,:)�� ,A(:,j)�з�����
%A(i,:)��ȡA�ĵ�i��
%A(:,j)��ȡA�ĵ�j��
clear;                                      %����ϴ������ɵ��ڴ��еĶ���
clc;                                        %����
m=input('please input  m=');                %����m��ֵ=8,16,32,64,128,256
alpha=input('please input  alpha=');        %����alpha��ֵ
h=1/(m+1);                                      %���񲽳�
n=m*m;                                      %���ɵľ���ά��
omega = pi;
mu = 1;
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
Umin = 1.00007
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if(alpha==0&&m==512)
eigS = eigs(T,W);
Umax = max(abs(eigS))
Umin = 1.00002
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

disp('LU�ֽ⣺')
[IT,res]=LU(W,T,n,alpha,p,q)


