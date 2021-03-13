function [IT,res]=cpcholesky(W,T,n,alpha,p,q)
r0=norm([p;q]);
gamma=(1.0/alpha);
P=symamd(W);
U=chol(W(P,P));
L=U';

res=10;
IT=0;
x=sparse(n,1);
y=sparse(n,1);
while(res>1e-6)&&(IT<200)
   t=T(P,P)*y+p(P);
   x=zuigan(L,t,n);
   t1=(1-gamma)*W(P,P)*y-gamma*T(P,P)*x+gamma*q(P);
   y=zuigan(L,t1,n);

   err1=p(P)-W(P,P)*x+T(P,P)*y;
   err2=q(P)-T(P,P)*x-W(P,P)*y;

   res=norm([err1;err2])/r0;
   IT=IT+1;
end
IT;
res;
end

function Z=zuigan(L,temp,n)
Z=sparse(n,1);
R=sparse(n,1);
R=L\temp;
Z=(L')\R;
end
