clc;
clear;
umin = 0.5;
umax = 3.5;
alphastep = 0.1;
for i=1:100
    alpha(i) = 0.8+i*alphastep;    
f(i)= 1-(1+umin*umin)/(alpha(i));
g(i) = 1-(1+umax*umax)/(alpha(i));
af(i)= abs(1-(1+umin*umin)/(alpha(i)));
ag(i) = abs(1-(1+umax*umax)/(alpha(i)));
end
figure(1)
plot(alpha,f,'r*')
hold on
plot(alpha,g,'b--')
hold on
xlabel('\alpha','FontSize',20);
ylabel('function value','FontSize',20);
legend({'f(\alpha)','g(\alpha)'},'Location','southeast','FontSize',20)

figure(2)
plot(alpha,af,'r*')
hold on
plot(alpha,ag,'b--')
hold on
xlabel('\alpha','FontSize',20);
ylabel('function absolute value','FontSize',20);
legend({'|f(\alpha)|','|g(\alpha)|'},'Location','northeast','FontSize',20)
