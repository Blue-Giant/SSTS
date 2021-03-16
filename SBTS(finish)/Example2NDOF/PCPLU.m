function [IT,res]=PCPLU(W,T,n,alpha,p,q)
r0=norm([p;q]);
P=symamd(W);
gamma = 1.0/alpha;

[L,U]=lu(W(P,P));

res=10;
IT=0;
x=sparse(n,1);
y=sparse(n,1);
   while(res>1e-6)&&(IT<500)
       t=T(P,P)*y+p(P);
       x=zuigan(L,U,t,n);
       t1=(1-gamma)*W(P,P)*y-gamma*T(P,P)*x+gamma*q(P);
       y=zuigan(L,U,t1,n);
       
       s=(1-gamma)*W(P,P)*y-gamma*T(P,P)*x+gamma*q(P);
%        s=(1-1/alpha)*((1-1/alpha)*W*y-(1/alpha)*T*x+(1/alpha)*q)-(1/alpha)*T*x+(1/alpha)*q;
       y=zuigan(L,U,s,n);
       s1=T(P,P)*y+p(P);
       x=zuigan(L,U,s1,n);
       
       err1=p(P)-W(P,P)*x+T(P,P)*y;
       err2=q(P)-T(P,P)*x-W(P,P)*y;
       
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

   