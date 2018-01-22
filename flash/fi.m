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
