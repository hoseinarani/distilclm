function [B,x,y,zx,zy]=twophaseLL(T,P,comp,K,n,k)
countrtz=50;
countrtbisec=100;
R=8.314;
er1=1;
er2=0;
a=1;
j=1;
kk=10;
while abs(er1-er2) > 1e-2

    for B1=eps:1/kk:1-eps
        if (RR(B1,k,n)*RR(B1+1/kk,k,n))<0
            B=B1+.5e-3;
        end
    end
    
    x=n./(B+k*(1-B));
    y=k.*x;
    zy=Zl(T,P,comp,K,y);
    zx=Zl(T,P,comp,K,x);
    for i=1:size(n,2)
        k(i)=fi(T,P,comp,K,x,zx,i)./fi(T,P,comp,K,y,zy,i);
    end
    er1=er2;
    er2=(B);
    j=j+1;
end
er1=1;
er2=0;
while abs(er1-er2) > 1e-5
    for i=1:size(n,2)
        e=[(k(i))-(fi(T,P,comp,K,x,zx,i))/(fi(T,P,comp,K,y,zy,i))
            RR(B,k,n)];
        d=1e-4;
        dRRdk=(RR(B,k+d,n)-RR(B,k,n))/d;
        dRRdB=(RR(B+d,k,n)-RR(B,k,n))/d;

        J=[  1         0
            dRRdk dRRdB ];
        delta=solution2(J,e);
        k(i)=(k(i))+delta(1);
        B=B+delta(2);
    end
    x=n./(B+k*(1-B));
    y=k.*x;
    zy=Zl(T,P,comp,K,y);
    zx=Zl(T,P,comp,K,x);
    for i=1:size(n,2)
        k(i)=fi(T,P,comp,K,x,zx,i)./fi(T,P,comp,K,y,zy,i);
    end
    er1=er2;
    er2=(norm(x)+norm(y)+B);
end

