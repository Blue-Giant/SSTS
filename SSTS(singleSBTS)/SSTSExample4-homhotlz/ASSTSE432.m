%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is  
%                 (W+iT)x = b
%     W =K+[(3-sqrt(3))/tau]I         T =K+[(3+sqrt(3))/tau]I
%          I(��λ��)                 K=I\kron Vm + Vm\kron I
%     Vm = (h^-2)*[tridiag(-1,2,-1)]   h=1/m        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           m=32ʱ���������������Ų����ı仯���������ĸı䣬��ͼ˵��          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��һ��6x6�����У���ȡ���ĵ�һ��Ԫ���γ�һ��6ά��������A(i,:)�� ,A(:,j)�з�����
%A(i,:)��ȡA�ĵ�i��
%A(:,j)��ȡA�ĵ�j��
clear;                                      %����ϴ������ɵ��ڴ��еĶ���
clc;                                        %����
m=32;                %����m��ֵ=8,16,32,64,128,256
sigma1=50;
sigma2=100;
h=1/(m+1);                                      %���񲽳�
n=m*m;                                      %���ɾ���ά��
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%���ԽǾ���(-1,2-1)
Vm=(1/(h*h))*V;                             
K=kron(speye(m),Vm)+kron(Vm,speye(m));


W=K+sigma1*speye(n);
W=h*h*W;
T=sigma2*speye(n);
T=h*h*T;
p=(W-T)*ones(n,1);
q=(W+T)*ones(n,1);



WW = W+T;
TT = T-W;
pp= p+q;
qq= q-p;

for i=1:20
   alpha(i) = 0.95+i*0.05;
   [IT(i),res(i)]=LU(WW,TT,n,alpha(i),pp,qq);
end

% plot(alpha,IT,'*r')
% xlabel('\alpha')
% ylabel('IT')
% title('ASSTS')
% set(get(gca,'title'),'FontSize',15,'FontName','Times New Roman');%���ñ��������С������  
% set(gca,'FontName','Times New Roman','FontSize',10)%���������������С������  
% set(get(gca,'XLabel'),'FontSize',15,'FontName','Times New Roman');%����X������������С������  
% set(get(gca,'YLabel'),'FontSize',15,'FontName','Times New Roman');%����Y������������С������ 
% axis([0.9 2 0 30])

plot(alpha,IT,'*b')
xlabel('\alpha')
ylabel('IT')
title('ASSTS')
set(get(gca,'title'),'FontSize',20,'FontName','Times New Roman','Color','k');%���ñ��������С������  
set(gca,'FontName','Times New Roman','FontSize',20)%���������������С������  
set(get(gca,'XLabel'),'FontSize',20,'FontName','Times New Roman','Color','k');%����X������������С������  
set(get(gca,'YLabel'),'FontSize',20,'FontName','Times New Roman','Color','k');%����Y������������С������ 
axis([0.9 2 0 100])


