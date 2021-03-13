function [IT,res]=PCPLU(W,T,n,alpha,p,q)
gamma=(1.0/alpha);
r0=norm([p;q]);
P=symamd(W);

[LP,UP]=lu(W(P,P));

res=10;
IT=0;
x=sparse(n,1);
y=sparse(n,1);
   while(res>1e-6)&&(IT<500)
       t=T(P,P)*y+p(P);
       x=zuigan(LP,UP,t,n);
       t1=(1.0-gamma)*W(P,P)*y-gamma*T(P,P)*x+gamma*q(P);
       y=zuigan(LP,UP,t1,n);
             
       err1=p(P)-W(P,P)*x+T(P,P)*y;
       err2=q(P)-T(P,P)*x-W(P,P)*y;
       
       res=norm([err1;err2])/r0;
       IT=IT+1;
   end
end

function Z=zuigan(L,U,temp,n)
Z=sparse(n,1);
R=sparse(n,1);
R=L\temp;
Z=U\R;
end

   