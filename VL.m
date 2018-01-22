function [B,x,y,zx,zy]=VL(T,P,z,comp)
K=zeros(size(z,2),size(z,2));
K=zeros(size(z,2),size(z,2));
R=8.314;
kz=zeros(size(z))+1e+6;
kz=(comp.Pc./P.*exp(5.37.*(1+comp.w).*(1-comp.Tc./T)));
[B,x,y,zx,zy]=twophaseVL(T,P,comp,K,z,kz)

function [B,x,y,zx,zy]=twophaseVL(T,P,comp,K,n,k)
R=8.314;
n1=.5;
er1=1;
er2=0;
a=1;
kk=50;
nextiter=1;
B=0;
if sum(n.*k)<=1
    VorL='L';
    %     nextiter=0;
elseif sum(n./k)<=1
    VorL='V';
    %     nextiter=0;
end
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
    x=x./norm(x);
    y=y./norm(y);
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

function zi=Zl(T,P,comp,K,y)
i=1;
R=8.314;
b=.07779*R.*comp.Tc./comp.Pc;
bm=sum(b.*y);
a=.45724.*(1+(.37464+1.54226.*comp.w-.26992.*comp.w.^2).*(1-(T./comp.Tc).^.5)).^2.*R^2.*comp.Tc.^2./comp.Pc;
am=0;
for i=1:size(y,2)
    for j=1:size(y,2)
        am=am+(a(j)*a(i))^.5*y(j)*y(i)*(1-K(i,j));
    end
end
B=bm*P/R/T;
q=am/bm/R/T;
a2=-(1-B);
a1=-(-q*B+2*B*(1+B)+B^2);
a0=-B^2*q+(1+B)*B^2;
% disp('no1')
% min(real(roots([1 a2 a1 a0])))
q=1/3*a1-1/9*a2^2;
r=1/6*(a1*a2-3*a0)-1/27*a2^3;
if q^3+r^2>0
    if r+(q^3+r^2)^.5<0
        s1=-abs(r+(q^3+r^2)^.5)^(1/3);
    else
        s1=(r+(q^3+r^2)^.5)^(1/3);
    end
    if r-(q^3+r^2)^.5<0
        s2=-abs(r-(q^3+r^2)^.5)^(1/3);
    else
        s2=(r-(q^3+r^2)^.5)^(1/3);
    end
    zi=[(s1+s2)-a2/3];
elseif q^3+r^2==0
    if r<0
        s1=-abs(r)^(1/3);
    else
        s1=(r)^(1/3);
    end
    if r<0
        s2=-abs(r)^(1/3);
    else
        s2=(r)^(1/3);
    end
    z=[(s1+s2)-a2/3
        -.5*(s1+s2)-a2/3];
    z=sort(z);
    zi=z(1);
    if z(1)<=0
        zi=z(2);
    end
else
    d=(-(q^3+r^2))^.5;
    b=(r^2+d^2)^.5;
    z=[2*b^(1/3)*cos(acos(r/b)/3)-a2/3
        -b^(1/3)*cos(acos(r/b)/3)-a2/3-3^.5*b^(1/3)*sin(acos(r/b)/3)
        -b^(1/3)*cos(acos(r/b)/3)-a2/3+3^.5*b^(1/3)*sin(acos(r/b)/3) ];
    z=sort(z);
    zi=z(1);
    if z(1)<=0
        zi=z(2);
    end
    if z(2)<=0
        zi=z(3);
    end
end
% disp('no2')
% zi
function zi=Zv(T,P,comp,K,y)
% global countrtz
% countrtz=100;
R=8.314;
b=.07779*R.*comp.Tc./comp.Pc;
bm=sum(b.*y);
a=.45724.*(1+(.37464+1.54226.*comp.w-.26992.*comp.w.^2).*(1-(T./comp.Tc).^.5)).^2.*R^2.*comp.Tc.^2./comp.Pc;
am=0;
for i=1:size(y,2)
    for j=1:size(y,2)
        am=am+(a(j)*a(i))^.5*y(j)*y(i)*(1-K(i,j));
    end
end
B=bm*P/R/T;
q=am/bm/R/T;
a2=-(1-B);
a1=-(-q*B+2*B*(1+B)+B^2);
a0=-B^2*q+(1+B)*B^2;
%disp('no1')
%((roots([1 a2 a1 a0])))
q=1/3*a1-1/9*a2^2;
r=1/6*(a1*a2-3*a0)-1/27*a2^3;

if q^3+r^2>0
    if r+(q^3+r^2)^.5<0
        s1=-abs(r+(q^3+r^2)^.5)^(1/3);
    else
        s1=(r+(q^3+r^2)^.5)^(1/3);
    end
    if r-(q^3+r^2)^.5<0
        s2=-abs(r-(q^3+r^2)^.5)^(1/3);
    else
        s2=(r-(q^3+r^2)^.5)^(1/3);
    end
    zi=[(s1+s2)-a2/3];
elseif q^3+r^2==0
    if r<0
        s1=-abs(r)^(1/3);
    else
        s1=(r)^(1/3);
    end
    if r<0
        s2=-abs(r)^(1/3);
    else
        s2=(r)^(1/3);
    end
    z=[(s1+s2)-a2/3
        -.5*(s1+s2)-a2/3];
    zi=max(z);
    
else
    d=(-(q^3+r^2))^.5;
    b=(r^2+d^2)^.5;
    z=[2*b^(1/3)*cos(acos(r/b)/3)-a2/3
        -b^(1/3)*cos(acos(r/b)/3)-a2/3-3^.5*b^(1/3)*sin(acos(r/b)/3)
        -b^(1/3)*cos(acos(r/b)/3)-a2/3+3^.5*b^(1/3)*sin(acos(r/b)/3) ];
    zi=max(z);
end
% disp('no2')
% zi


function ansfi=fi(T,P,comp,K,y,z,n)
R=8.314;
b=.07779*R*comp.Tc./comp.Pc;
bm=sum(b.*y);
a=.45724.*(1+(.37464+1.54226*comp.w-.26992*comp.w.^2).*(1-(T./comp.Tc).^.5)).^2.*R^2.*comp.Tc.^2./comp.Pc;
am=0;
for i=1:size(y,2)
    for j=1:size(y,2)
        am=am+(a(j)*a(i))^.5*y(j)*y(i)*(1-K(i,j));
    end
end
B=bm*P/R/T;
q=am/bm/R/T;
ai=0;
for j=1:size(y,2)
    ai=ai+y(j)*(a(n)*a(j))^.5*(1-K(n,j));
end
% if (z-B)<eps
%   10
% end
ai=2*ai-am;
qi=q*(1+ai/am-b(n)/bm);
I=1/2/2^.5*log(abs((z+(1+2^.5)*B)/(z+(1-2^.5)*B)));
ansfi=exp(b(n)/bm*(z-1)-log(abs(z-B))-qi*I);

function x=solution2(J,e)
eror1=1;
eror2=0;
x=[0 0];
while abs(eror1-eror2)>1e-6
    eror2=eror1;
    x(1)=(-e(1)-(J(1,2)*x(2)))/J(1,1);
    x(2)=(-e(2)-(J(2,1)*x(1)))/J(2,2);
    eror1=norm(x);
end

function R=RR(B,k,n)
R=sum(n.*(k-1)./(B+k.*(1-B)));
