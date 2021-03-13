function [IT,res]=LU(WW,TT,n,alpha,p,q)
r0=norm([p;q]);

[L,U]=lu(WW);

res=10;
IT=0;
x=sparse(n,1);
y=sparse(n,1);
% tic;
   while(res>1e-6)&&(IT<200)
       s=TT*y+p;
       x=zuigan(L,U,s,n);
       t=(1-1/alpha)*WW*y-(1/alpha)*TT*x+(1/alpha)*q;
       y=zuigan(L,U,t,n);       
       err1=p-WW*x+TT*y;
       err2=q-TT*x-WW*y;       
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

   