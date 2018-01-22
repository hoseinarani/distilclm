function ps=psat(T,P,Tc,Pc,w,K)
global countrt ps
countrt=500;
dp=1e-4;
ps=ones(size(w)).*P;
for i=1:size(w,2)
    fl=0;
    fv=1;
    if T<Tc(i)
        while abs(fl/fv-1)>1e-4
            %         ps1=ps(i);
            zv=Zv(T,ps(i),Tc(i),Pc(i),w(i),K,1);
            zl=Zl(T,ps(i),Tc(i),Pc(i),w(i),K,1);
            fl=fi(T,ps(i),Tc(i),Pc(i),w(i),K,1,zl,1);
            fv=fi(T,ps(i),Tc(i),Pc(i),w(i),K,1,zv,1);
            ps(i)=ps(i)*fl/fv;
%                     dfdp=(ff(T,ps(i)+dp,Tc(i),Pc(i),w(i),K,i)-ff(T,ps(i),Tc(i),Pc(i),w(i),K,i))/dp;
%                     ps(i)=ps(i)-ff(T,ps(i),Tc(i),Pc(i),w(i),K,i)/dfdp;

        end
    else
        ps(i)=Pc(i);
    end
end

function f=ff(T,P,Tc,Pc,w,K,i)
global ps
zv=Zv(T,P,Tc,Pc,w,K,1);
zl=Zl(T,P,Tc,Pc,w,K,1);
while zv==-1 | zl==-1
    if zv==-1
        ps(i)=ps(i)/2
        zv=Zv(T,P,Tc,Pc,w,K,1);
        zl=Zl(T,P,Tc,Pc,w,K,1);
    elseif zl==-1
        ps(i)=ps(i)*sqrt(2)
        zv=Zv(T,P,Tc,Pc,w,K,1);
        zl=Zl(T,P,Tc,Pc,w,K,1);
    end
end
f=fi(T,P,Tc,Pc,w,K,1,zv,1)-fi(T,P,Tc,Pc,w,K,1,zl,1);

