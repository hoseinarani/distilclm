
function [B,x,y,zx,zy]=twophaseVL(T,P,comp,K,n,k)
R=8.314;
n1=.5;
er1=1;
er2=0;
a=1;
kk=10;
nextiter=1;
B=0.5;
% if sum(n.*k)<=1
%     VorL='L';
%     nextiter=0;
% elseif sum(n./k)<=1
%     VorL='V';
%     nextiter=0;
% end
% if  nextiter==0
%     er1=er2;
% end
% % jj=0;
% % countrtbisec=50;
while abs(er1-er2) > 1e-3
%         jj=jj+1;
%     if jj>countrtbisec
%         ans3=0;
%         return
%     end
    for B1=eps:1/kk:1-eps
        if (RR(B1,k,n)*RR(B1+1/kk,k,n))<0
            B=B1+.5e-3;
        end
    end
    x=n./(B+k*(1-B));
    y=k.*x;
    zy=Zv(T,P,comp,K,y);
    zx=Zl(T,P,comp,K,x);
    for i=1:size(n,2)
        k(i)=fi(T,P,comp,K,x,zx,i)./fi(T,P,comp,K,y,zy,i);
    end
    er1=er2;
    er2=(B);
    if sum(n.*k)<=1 | B>=1 %& sum(n./k)<=1) | B<=0
        VorL='L';
        er1=er2;
        nextiter=0;
    elseif sum(n./k)<=1 | B<=0
        VorL='V';
        er1=er2;
        nextiter=0;
    end
end
if  nextiter==0
    er1=er2;
else
    er1=1;
    er2=0;
end
count=1;
while abs(er1-er2) > 1e-5
%     if jj>counter | B<0 |  x<0 | y<0 | 1-(B)<0 | B>1 | x>1 | y>1  | 1-(B)>1
%         ans2=0;
%         return
%     end
    for i=1:size(n,2)
        e=[(k(i))-(fi(T,P,comp,K,x,zx,i))/(fi(T,P,comp,K,y,zy,i))
            RR(B,k,n)];
        d=zeros(1,size(n,2));
        d(1,i)=1e-4;
        dRRdk=(RR(B,k+d,n)-RR(B,k,n))/d(i);
        dRRdB=(RR(B+d(i),k,n)-RR(B,k,n))/d(i);

        J=[  1         0
            dRRdk  dRRdB ];

        delta=solution2(J,e);
        k(i)=(k(i))+delta(1);
        B=B+delta(2);
    end
    x=n./(B+k*(1-B));
    y=k.*x;
    zy=Zv(T,P,comp,K,y);
    zx=Zl(T,P,comp,K,x);
    for i=1:size(n,2)
        k(i)=fi(T,P,comp,K,x,zx,i)./fi(T,P,comp,K,y,zy,i);
    end
    er1=er2;
    er2=(norm(x)+norm(y)+B);
    VorL='VL';
    if sum(n.*k)<=1 | B>=1 %& sum(n./k)<=1) | B<=0
        VorL='L';
        er1=er2;
    elseif sum(n./k)<=1 | B<=0
        VorL='V';
        er1=er2;
    end
    count=count+1;
end
if VorL=='L'
    B=1;
    x=n;
    y=zeros(1,size(n,2));
    zx=Zl(T,P,comp,K,x);
    zy=0;
elseif VorL=='V'
    B=0;
    x=zeros(1,size(n,2));
    y=n;
    zx=0;
    zy=Zv(T,P,comp,K,y);
end






