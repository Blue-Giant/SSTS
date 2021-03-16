function [IT,res]=LU(W,T,n,alpha)
p=(W-T)*ones(n,1);
q=(W+T)*ones(n,1);
r0=norm([p,q]);

[L,U]=lu(W);

res=10;
IT=0;
x=sparse(n,1);
y=sparse(n,1);
tic;
   while(res>1e-6)&&(IT<1000)
       t=T*y+p;
       x=zuigan(L,U,t,n);
%        t1=(1-1/alpha)*W*y-(1/alpha)*T*x+(1/alpha)*q;
%        y=zuigan(L,U,t1,n);
%        s=(1-1/alpha)*W*y-(1/alpha)*T*x+(1/alpha)*q;
       s=(1-1/alpha)*((1-1/alpha)*W*y-(1/alpha)*T*x+(1/alpha)*q)-(1/alpha)*T*x+(1/alpha)*q;
       y=zuigan(L,U,s,n);
       s1=T*y+p;
       x=zuigan(L,U,s1,n);
       
       err1=p-W*x+T*y;
       err2=q-T*x-W*y;
       
       res=norm([err1,err2])/r0;
       IT=IT+1;
   end
   toc;
IT;
res;
end

function Z=zuigan(L,U,temp,n)
Z=sparse(n,1);
R=sparse(n,1);
R=L\temp;
Z=U\R;
end

   