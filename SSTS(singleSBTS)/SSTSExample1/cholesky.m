function [IT,res]=cholesky(W,T,n,alpha,p, q)
r0=norm([p,q]);
gamma = 1.0/alpha;
U=chol(W);
L=U';

res=10;
IT=0;
x=sparse(n,1);
y=sparse(n,1);
   while(res>1e-6)&&(IT<200)
       t=T*y+p;
       x=zuigan(L,t,n);
       t1=(1-gamma)*W*y-gamma*T*x+gamma*q;
       y=zuigan(L,t1,n);       
       err1=p-W*x+T*y;
       err2=q-T*x-W*y;
       
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
