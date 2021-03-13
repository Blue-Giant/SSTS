function [IT,res]=cpcholesky(W,T,n,alpha,b)
P=symamd(W);
r0=norm(b);
I=speye(n);

Uw=chol(alpha*I+W(P,P));
Lw=Uw';
Ut=chol(alpha*I+T(P,P));
Lt=Ut';

res=10;
IT=0;
u=sparse(n,1);
% tic;
   while(res>1e-6)&&(IT<1000)
       t=(alpha*I-1i*T(P,P))*u+b(P);
       u=zuigan(Lw,t,n);
       t1=(alpha*I+1i*W(P,P))*u-1i*b(P);
       u=zuigan(Lt,t1,n);
       error=b(P)-(W(P,P)+1i*T(P,P))*u;       
       res=norm(error)/r0;
       IT=IT+1;
   end
%    toc;
IT;
res;
end

function Z=zuigan(L,temp,n)
Z=sparse(n,1);
R=sparse(n,1);
R=L\temp;
Z=(L')\R;
end
