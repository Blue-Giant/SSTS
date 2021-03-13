function [IT,res]=PCPLU(W,T,n,alpha,b)
r0=norm(b);
P=symamd(W);
I=speye(n);

[Lw,Uw]=lu(alpha*I+W(P,P));
[Lt,Ut]=lu(alpha*I+T(P,P));

res=10;
IT=0;
u=sparse(n,1);
   while(res>1e-6)&&(IT<1000)
       t=(alpha*I-1i*T(P,P))*u+b(P);
       u=zuigan(Lw,Uw,t,n);
       t1=(alpha*I+1i*W(P,P))*u-1i*b(P);
       u=zuigan(Lt,Ut,t1,n);
       error=b(P)-(W(P,P)+1i*T(P,P))*u;       
       res=norm(error)/r0;
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

   