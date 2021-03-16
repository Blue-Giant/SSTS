function [IT,res]=LU(WW,TT,n,alpha,pp,qq)
r0=norm([pp;qq]);

[L,U]=lu(WW);

res=10;
IT=0;
x = sparse(n,1);
y = sparse(n,1);
   while(res>1e-6)&&(IT<200)
       t=TT*y+pp;
       x = zuigan(L,U,t,n);
       t1 =( 1-1/alpha)*WW*y-(1/alpha)*TT*x+(1/alpha)*qq;
       y=zuigan(L,U,t1,n);
       
       s=(1-1/alpha)*WW*y-(1/alpha)*TT*x+(1/alpha)*qq;
       y=zuigan(L,U,s,n);
       s1=TT*y+pp;
       x=zuigan(L,U,s1,n);
       
       err1=pp-WW*x+TT*y;
       err2=qq-TT*x-WW*y;
       
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

   