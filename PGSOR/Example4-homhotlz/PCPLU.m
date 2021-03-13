function [IT,res]=PCPLU(W,T,n,alpha,p,q)
r0=norm([p;q]);
P  = symamd(W);
[L,U]=lu(W(P,P));
res=10;
IT=0;
x=sparse(n,1);
y=sparse(n,1);
% tic;
   while(res>1e-6)&&(IT<200)
       t=(1-alpha)*W(P,P)*x  +alpha*T(P,P)*y + alpha*p(P);
       x=zuigan(L,U,t,n);       
       s=(1-alpha)*W(P,P)*y - (alpha)*T(P,P)*x + (alpha)*q(P);
       y=zuigan(L,U,s,n);       
       err1=p(P)-W(P,P)*x+T(P,P)*y;
       err2=q(P)-T(P,P)*x-W(P,P)*y;       
       res=norm([err1;err2])/r0;
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

   