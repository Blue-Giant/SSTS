%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The complex symmetric linear systems of equations is  
%                 (W+iT)x = b
%     W =K+[(3-sqrt(3))/tau]I         T =K+[(3+sqrt(3))/tau]I
%          I(单位阵)                 K=I\kron Vm + Vm\kron I
%     Vm = (h^-2)*[tridiag(-1,2,-1)]   h=1/m        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           m=32时，迭代步数会随着参数的变化发生怎样的改变，画图说明          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%从一个6x6矩阵中，提取它的第一行元素形成一个6维行向量。A(i,:)行 ,A(:,j)列方法。
%A(i,:)提取A的第i行
%A(:,j)提取A的第j列
clear;                                      %清除上次所生成的内存中的东西
clc;                                        %清屏
m=32;                %输入m的值=8,16,32,64,128,256
sigma1=50;
sigma2=100;
h=1/(m+1);                                      %网格步长
n=m*m;                                      %生成矩阵维数
V=diag(-1*ones(m-1,1),1)+2*speye(m)+diag(-1*ones(m-1,1),-1);%三对角矩阵(-1,2-1)
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
% set(get(gca,'title'),'FontSize',15,'FontName','Times New Roman');%设置标题字体大小，字型  
% set(gca,'FontName','Times New Roman','FontSize',10)%设置坐标轴字体大小，字型  
% set(get(gca,'XLabel'),'FontSize',15,'FontName','Times New Roman');%设置X坐标标题字体大小，字型  
% set(get(gca,'YLabel'),'FontSize',15,'FontName','Times New Roman');%设置Y坐标标题字体大小，字型 
% axis([0.9 2 0 30])

plot(alpha,IT,'*b')
xlabel('\alpha')
ylabel('IT')
title('ASSTS')
set(get(gca,'title'),'FontSize',20,'FontName','Times New Roman','Color','k');%设置标题字体大小，字型  
set(gca,'FontName','Times New Roman','FontSize',20)%设置坐标轴字体大小，字型  
set(get(gca,'XLabel'),'FontSize',20,'FontName','Times New Roman','Color','k');%设置X坐标标题字体大小，字型  
set(get(gca,'YLabel'),'FontSize',20,'FontName','Times New Roman','Color','k');%设置Y坐标标题字体大小，字型 
axis([0.9 2 0 100])


