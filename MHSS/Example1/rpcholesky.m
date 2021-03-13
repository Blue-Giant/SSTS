function [IT,res]=rpcholesky(W,T,n,alpha,b)
R=symamd(W);
r0=norm(b);
I=speye(n);

Uw=chol(alpha*I+W(R,R));
Lw=Uw';
Ut=chol(alpha*I+T(R,R));
Lt=Ut';

res=10;
IT=0;
u=sparse(n,1);
tic;
   while(res>1e-6)&&(IT<1000)
       t=(alpha*I-1i*T(R,R))*u+b(R);
       u=zuigan(Lw,t,n);
       t1=(alpha*I+1i*W(R,R))*u-1i*b(R);
       u=zuigan(Lt,t1,n);
       error=b(R)-(W(R,R)+1i*T(R,R))*u;       
       res=norm(error)/r0;
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
