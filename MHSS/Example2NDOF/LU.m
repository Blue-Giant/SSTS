function [IT,res]=LU(W,T,n,alpha,b)
r0=norm(b);
I=speye(n);

[Lw,Uw]=lu(alpha*I+W);
[Lt,Ut]=lu(alpha*I+T);

res=10;
IT=0;
u=sparse(n,1);
% tic;
   while(res>1e-6)&&(IT<1000)
       t=(alpha*I-1i*T)*u+b;
       u=zuigan(Lw,Uw,t,n);
       t1=(alpha*I+1i*W)*u-1i*b;
       u=zuigan(Lt,Ut,t1,n);
       error=b-(W+1i*T)*u;       
       res=norm(error)/r0;
       IT=IT+1;
   end
%    toc;
IT;
res;
end

function Z=zuigan(L,U,temp,n)
Z=sparse(n,1);
R=sparse(n,1);
R=L\temp;
Z=U\R;
end

   