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
