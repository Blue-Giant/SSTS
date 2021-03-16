function [IT,res]=cpcholesky(W,T,n,alpha)
P=symamd(W);
p=(W(P,P)-T(P,P))*ones(n,1);
q=(W(P,P)+T(P,P))*ones(n,1);
r0=norm([p,q]);

U=chol(W(P,P));
L=U';

res=10;
IT=0;
x=sparse(n,1);
y=sparse(n,1);
tic;
   while(res>1e-6)&&(IT<200)
       t=T(P,P)*y+p;
       x=zuigan(L,t,n);
       t1=(1-1/alpha)*W(P,P)*y-(1/alpha)*T(P,P)*x+(1/alpha)*q;
       y=zuigan(L,t1,n);
       
       s=(1-1/alpha)*W(P,P)*y-(1/alpha)*T(P,P)*x+(1/alpha)*q;
       y=zuigan(L,s,n);
       s1=T(P,P)*y+p;
       x=zuigan(L,s1,n);
       
       err1=p-W(P,P)*x+T(P,P)*y;
       err2=q-T(P,P)*x-W(P,P)*y;
       
       res=norm([err1,err2])/r0;
       IT=IT+1;
   end
   toc;
IT;
res;
end

function Z=zuigan(L,temp,n)
Z=sparse(n,1);
R=sparse(n,1);
R=L\temp;
Z=(L')\R;
end
