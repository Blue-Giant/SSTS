function [IT,res]=LU(W,T,n,alpha,p,q)
r0=norm([p;q]);
[L,U]=lu(W);
res=10;
IT=0;
x=sparse(n,1);
y=sparse(n,1);
   while(res>1e-6)&&(IT<200)
       t=(1-alpha)*W*x  +alpha*T*y + alpha*p;
       x=zuigan(L,U,t,n);       
       s=(1-alpha)*W*y - (alpha)*T*x + (alpha)*q;
       y=zuigan(L,U,s,n);       
       err1=p-W*x+T*y;
       err2=q-T*x-W*y;       
       res=norm([err1;err2])/r0;
       IT=IT+1;
   end
IT;
res;
end

function Z=zuigan(L,U,temp,n)
Z=sparse(n,1);
R=sparse(n,1);
R=L\temp;
Z=U\R;
end

   