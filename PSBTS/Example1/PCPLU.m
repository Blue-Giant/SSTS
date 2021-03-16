function [IT,res]=PCPLU(WW,TT,n,alpha,pp,qq)
r0=norm([pp;qq]);
P=symamd(WW);
[L,U]=lu(WW(P,P));

res=10;
IT=0;
x = sparse(n,1);
y = sparse(n,1);
while(res>1e-6)&&(IT<200)
   t=TT(P,P)*y+pp(P);
   x = zuigan(L,U,t,n);
   t1 =( 1-1/alpha)*WW(P,P)*y-(1/alpha)*TT(P,P)*x+(1/alpha)*qq(P);
   y=zuigan(L,U,t1,n);

   s=(1-1/alpha)*WW(P,P)*y-(1/alpha)*TT(P,P)*x+(1/alpha)*qq(P);
   y=zuigan(L,U,s,n);
   s1=TT(P,P)*y+pp(P);
   x=zuigan(L,U,s1,n);

   err1=pp(P)-WW(P,P)*x+TT(P,P)*y;
   err2=qq(P)-TT(P,P)*x-WW(P,P)*y;

   res=norm([err1;err2])/r0;
   IT=IT+1;
end
end

function Z=zuigan(L,U,temp,n)
Z=sparse(n,1);
R=sparse(n,1);
R=L\temp;
Z=U\R;
end

   