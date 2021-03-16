function [IT,res]=zxdcholesky(W,T,n,alpha)
D = symmmd(W);
p = (W(D,D)-T(D,D))*ones(n,1);
q = (W(D,D)+T(D,D))*ones(n,1);
r0 = norm([p,q]);

U=chol(W(D,D));
L=U';

res=10;
IT=0;
x=sparse(n,1);
y=sparse(n,1);
tic;
   while(res>1e-6)&&(IT<200)
       t=T(D,D)*y+p;
       x=zuigan(L,t,n);
       t1=(1-1/alpha)*W(D,D)*y-(1/alpha)*T(D,D)*x+(1/alpha)*q;
       y=zuigan(L,t1,n);
       
       s=(1-1/alpha)*W(D,D)*y-(1/alpha)*T(D,D)*x+(1/alpha)*q;
       y=zuigan(L,s,n);
       s1=T(D,D)*y+p;
       x=zuigan(L,s1,n);
       
       err1=p-W(D,D)*x+T(D,D)*y;
       err2=q-T(D,D)*x-W(D,D)*y;
       
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
