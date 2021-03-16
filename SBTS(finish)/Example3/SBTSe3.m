%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is          
%                           (W+iT)*u = b; 
%     W=10(I\kron Vc + Vc\kron I)+9(e1*em'+em*e1')     T=I\kron Vc + Vc\kron I
%     V = tridiag(-1,2,-1)   Vc = V -  e1*em'+em*e1'   
%     e1 ��λ����(1,0,0,0,......,,0)    em ��λ����(0,0,0,0,......,,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SBTS�����Ĳ���ѡȡ�б�
%m            8          16         32         64         128        256  
%alpha      1.3747     1.7470     2.8821     6.8786     22.1690     82.3060
%beta       0.7858     0.7005     0.6050     0.5392     0.5115      0.5031   
%Umax       0.3962     0.6667     1.2183     2.3270     4.5473      8.9892
%Umin       0.0599     0.0551     0.0526     0.0513     0.0507      0.0500
%IT            5         8          16         43        146         560
%RES        9.5624     5.6641     6.5896     8.4509     9.5010      9.7746
%-----------------------------ȡ��λ��Ч����-------------------------------
%alpha(IT)  1.37(5)    1.75(8)    2.88(16)   6.88(43)   22.17(146)     48
%res        9.4695     5.8455     6.5230     8.4562     9.4762
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��һ��6x6�����У���ȡ���ĵ�һ��Ԫ���γ�һ��6ά��������A(i,:)�� ,A(:,j)�з�����
%A(i,:)��ȡA�ĵ�i��
%A(:,j)��ȡA�ĵ�j��
clear;                                      %����ϴ������ɵ��ڴ��еĶ���
clc;                                        %����
m=input('please input  m=');                %��8��m��ֵ=8,16,32,64,128,256
alpha=input('please input  alpha=');        %���������alpha��ֵ,���alpha=0,������Զ�����alpha��ֵ
h=1/m;                                      %���񲽳�
n=m*m;                                      %���ɾ���ά��
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%���ԽǾ���(1-,2-1)
e=speye(m);                                 %��ɢϡ��洢�ĵ�λ����
e1=e(:,1);                                  %��λ����(1,0,0,0,......,,0)
em=e(:,m);                                  %��λ����(0,0,0,0,......,,1)
Ve=e1*em'+em*e1';
Vc=V-Ve;

W=10*(kron(e,Vc)+kron(Vc,e)) +9*kron(Ve,e); %�Գ���������
W=h*h*W;                             
T=kron(e,V)+kron(V,e);                      %�Գƾ���
T=h*h*T;

%����ֿ�����mat2cell
% mat2cell(W,())
%�����ֵ
if (alpha==0&&m~=256)
S=W\T;
eigS = eig(full(S));
Umax = max(abs(eigS))
Umin = min(abs(eigS))
temp = 2+ Umax*Umax + Umin*Umin;
alpha = (temp+sqrt(temp*temp-temp-temp))/2
beta = (temp-sqrt(temp*temp-temp-temp))/2
end

if (alpha==0&&m==256)
eigS = eigs(T,W);
Umax = max(abs(eigS))
Umin = 0.05
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


