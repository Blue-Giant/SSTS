function [IT,res]=rpcholesky(W,T,n,alpha,p,q)
R=symamd(W);
r0=norm([p,q]);

U=chol(W(R,R));
L=U';

res=10;
IT=0;
x=sparse(n,1);
y=sparse(n,1);
   while(res>1e-6)&&(IT<200)
       t=T(R,R)*y+p(R);
       x=zuigan(L,t,n);
       t1=(1-1/alpha)*W(R,R)*y-(1/alpha)*T(R,R)*x+(1/alpha)*q(R);
       y=zuigan(L,t1,n);
       
       s=(1-1/alpha)*W(R,R)*y-(1/alpha)*T(R,R)*x+(1/alpha)*q(R);
       y=zuigan(L,s,n);
       s1=T(R,R)*y+p(R);
       x=zuigan(L,s1,n);
       
       err1=p(R)-W(R,R)*x+T(R,R)*y;
       err2=q(R)-T(R,R)*x-W(R,R)*y;
       
       res=norm([err1,err2])/r0;
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
