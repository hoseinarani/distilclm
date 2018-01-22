function [Lj,Vj,xji,yji,Tjk,Qj,Uj,Wj,Fj,Pj,N,nc]=columncalc(clm,comp,feed)
format short g;

% components C3, n-C4, n-C5
%consider that all feed in liquid
nc = size(comp.Mw,2);
Tc=comp.Tc;
Pc=comp.Pc;
w=comp.w;
K=zeros(nc);
N=clm.n_tray+2;%stages
NF=clm.f_tray+1;%feed stage
D=clm.distil_rate;
dp=clm.dp*1000; %pa/tray
rr=clm.rr;

TF=str2double(cell2mat(feed.data(1,2)))+273; %K
PF=str2double(cell2mat(feed.data(2,2)))*1000; %pa
P=str2double(cell2mat(feed.data(2,2)))*1000; %pa pressure of feed fo anthalpy calc
zFji=zeros(N,nc);
zFji(NF,:)=str2double(cellfun(@cellstr,feed.comp(:,2)))';
Fj=zeros(N,1);
Fj(NF)=str2double(cell2mat(feed.data(3,2)));

TFj=zeros(N,1)+TF;

%TFj(NF)=TF;
Uj=zeros(N,1);
Pj=zeros(N,1)+PF;%pa
for i=NF:N-1
    Pj(i+1)=Pj(i)+dp;
end
for i=NF:-1:2
    Pj(i-1)=Pj(i)-dp;
end
Uj(1)=D;
Wj=zeros(N,1);

Qj=zeros(N,1);


%initial guess
zF=zFji(NF,:);
Tjinit=TF;
Vjinit=Fj(NF)/2;
Ljinit=Fj(NF)/2;
Tjk=zeros(N,1)+Tjinit;%K
Vj=zeros(N,1)+Vjinit;
Lj=zeros(N,1)+Ljinit;
Vj(1)=0; % no vapor out of condenser
Tjk(:,2)=Tjk(:,1)+2;
% Tjk(1,:)=-17+273;
% Tjk(5,:)=93+273;
k=2;
for j=1:N
    xji(j,:)=zF';
    yji(j,:)=zF';
end

while abs(sum(Tjk(:,k)-Tjk(:,k-1)))>1
    abs(sum(Tjk(:,k)-Tjk(:,k-1)))
    pause(0.1)
    for j=1:N
        for i=1:nc
            zy=Zv(Tjk(j,k),Pj(j),comp,K,yji(j,i));
            zx=Zl(Tjk(j,k),Pj(j),comp,K,xji(j,i));
            Kji(j,i)=fi(Tjk(j,k),Pj(j),comp,K,xji(j,i),zx,i)/fi(Tjk(j,k),Pj(j),comp,K,yji(j,i),zy,i);
            if isnan(Kji(j,i))
                %hossein added: to ignore error,
                Kji(j,i)=1; 
            end
                
        end
    end
    for j=1:N
        
        if j>=2
            Aj(j)=Lj(j-1);% from algoritmo
%             Aj(j)=Vj(j)+sum(Fj(1:j-1)-Wj(1:j-1)-Uj(1:j-1))-Vj(1);
        end
        if j==N
            Bji(j,:)=-(Lj(j)+Uj(j)+(Vj(j)+Wj(j))*Kji(j,:));%% from algoritmo
%             Bji(j,:)=-(0+sum(Fj(1:j)-Wj(1:j)-Uj(1:j))-Vj(1)+Uj(j)+(Vj(j)+Wj(j))*Kji(j,:));
            Cji(j,:)=0;
        else
            Bji(j,:)=-(Lj(j)+Uj(j)+(Vj(j)+Wj(j))*Kji(j,:));%
%             Bji(j,:)=-(Vj(j+1)+sum(Fj(1:j)-Wj(1:j)-Uj(1:j))-Vj(1)+Uj(j)+(Vj(j)+Wj(j))*Kji(j,:));
            Cji(j,:)=Vj(j+1)*Kji(j+1,:);
        end
        
        Dji(j,:)=-Fj(j)*zFji(j,:);
        
    end
    for i=1:nc
        clear A B
        A=full(gallery('tridiag',Aj(2:end),Bji(:,i)',Cji(1:end-1,i)'));
        B=Dji(:,i);
        xji(:,i)=linsolve(A,B);
    end
    for j=1:N
        xji(j,:)=abs(xji(j,:))./sum(abs(xji(j,:)));
    end
%            xji(1,1)=.8;

    for j=1:N
        [Tjk(j,k+1),yji(j,:)]=bublepoint(Tjk(j,k),xji(j,:),Pj(j),comp,K);
    end
    Qj(1)=Vj(2)*Hv(comp,yji(2,:)',Tjk(2,k+1))+Fj(1)*Hl(comp,zFji(1,:)',TFj(1),P)-(Lj(1)+Uj(1))*Hl(comp,xji(1,:)',Tjk(1,k+1),P)...
        -(Vj(1)+Wj(1))*Hv(comp,yji(1,:)',Tjk(1,k+1));
    
    Qj(N)=(-Vj(1)*Hv(comp,yji(1,:)',Tjk(1,k+1))...
                -Lj(N)*Hl(comp,xji(N,:)',Tjk(N,k+1),P));
    for j=1:N
        if j==N
            Qj(N)=Qj(N)+Fj(j)*Hl(comp,zFji(j,:)',TFj(j),P)...
                -Uj(j)*Hl(comp,xji(j,:)',Tjk(j,k+1),P)...
                -Wj(j)*Hv(comp,yji(j,:)',Tjk(j,k+1));
        else
            Qj(N)=Qj(N)+Fj(j)*Hl(comp,zFji(j,:)',TFj(j),P)...
                -Uj(j)*Hl(comp,xji(j,:)',Tjk(j,k+1),P)...
                -Wj(j)*Hv(comp,yji(j,:)',Tjk(j,k+1))...
                -Qj(j);
        end
    end
    Qj(N)=(Qj(N));
    for j=2:N
        if j==2
            alfaj(j)=Hl(comp,xji(j-1,:)',Tjk(j-1,k+1),P)-Hv(comp,yji(j,:)',Tjk(j,k+1));
            betaj(j)=Hv(comp,yji(j+1,:)',Tjk(j+1,k+1))-Hl(comp,xji(j,:)',Tjk(j,k+1),P);
            gamaj(j)=(sum(Fj(1:j-1)-Wj(1:j-1)-Uj(1:j-1))-Vj(1))*(Hl(comp,xji(j,:)',Tjk(j,k+1),P)-Hl(comp,xji(j-1,:)',Tjk(j-1,k+1),P))...
                +Fj(j)*(Hl(comp,xji(j,:)',Tjk(j,k+1),P)-Hl(comp,zFji(j,:)',TFj(j),P))...
                +Wj(j)*(Hv(comp,yji(j,:)',Tjk(j,k+1))-Hl(comp,xji(j,:)',Tjk(j,k+1),P))+Qj(j);
             %Vj(j)=Lj(j)+Vj(1)-sum(Fj(1:j)-Wj(1:j)-Uj(1:j));
             Vj(j)=Lj(1)+Vj(1)+Uj(1)-Fj(1);
        elseif j==N
            alfaj(j)=Hl(comp,xji(j-1,:)',Tjk(j-1,k+1),P)-Hv(comp,yji(j,:)',Tjk(j,k+1));
            betaj(j)=-Hl(comp,xji(j,:)',Tjk(j,k+1),P);
            gamaj(j)=(sum(Fj(1:j-1)-Wj(1:j-1)-Uj(1:j-1))-Vj(1))*(Hl(comp,xji(j,:)',Tjk(j,k+1),P)-Hl(comp,xji(j-1,:)',Tjk(j-1,k+1),P))...
                +Fj(j)*(Hl(comp,xji(j,:)',Tjk(j,k+1),P)-Hl(comp,zFji(j,:)',TFj(j),P))...
                +Wj(j)*(Hv(comp,yji(j,:)',Tjk(j,k+1))-Hl(comp,xji(j,:)',Tjk(j,k+1),P))+Qj(j);
            Vj(j)=(gamaj(j-1)-alfaj(j-1)*Vj(j-1))/betaj(j-1);
        else
            alfaj(j)=Hl(comp,xji(j-1,:)',Tjk(j-1,k+1),P)-Hv(comp,yji(j,:)',Tjk(j,k+1));
            betaj(j)=Hv(comp,yji(j+1,:)',Tjk(j+1,k+1))-Hl(comp,xji(j,:)',Tjk(j,k+1),P);
            gamaj(j)=(sum(Fj(1:j-1)-Wj(1:j-1)-Uj(1:j-1))-Vj(1))*(Hl(comp,xji(j,:)',Tjk(j,k+1),P)-Hl(comp,xji(j-1,:)',Tjk(j-1,k+1),P))...
                +Fj(j)*(Hl(comp,xji(j,:)',Tjk(j,k+1),P)-Hl(comp,zFji(j,:)',TFj(j),P))...
                +Wj(j)*(Hv(comp,yji(j,:)',Tjk(j,k+1))-Hl(comp,xji(j,:)',Tjk(j,k+1),P))+Qj(j);
            Vj(j)=(gamaj(j-1)-alfaj(j-1)*Vj(j-1))/betaj(j-1);
        end
%     Vj=abs(Vj);% hossein added

    end
    for j=1:N
        if j==N
            Lj(j)=0+sum(Fj(1:j)-Wj(1:j)-Uj(1:j))-Vj(1);
        elseif j==1 % reflux ratio
            Lj(j)=rr*D;
%             Lj(j)=Vj(j+1)+sum(Fj(1:j)-Wj(1:j)-Uj(1:j))-Vj(1);

        else
            Lj(j)=Vj(j+1)+sum(Fj(1:j)-Wj(1:j)-Uj(1:j))-Vj(1);
        end;
    end
%     Lj=abs(Lj);% hossein added
    k=k+1;
end
Tj=Tjk(j,end);
clc







function H = Hm(coef,x,T)%J/mol
n=length(x);
for i=1:n
    Hi(i,:)=polyval(coef(i,:),T);
end;
H=sum(x.*Hi);
function [Tbp,y]=bublepoint(T,x,P,copm,K)
% options=optimset('Display','iter');options
answ=fsolve(@bublepointfunc,[T;x'],'',x,P,copm,K);
Tbp=answ(1);
y=answ(2:end);
function answ=bublepointfunc(var,x,P,comp,K)
T=var(1);
y=var(2:end)';

zy=Zv(T,P,comp,K,y);
zx=Zl(T,P,comp,K,x);
for i=1:size(x,2)
    ansf(i)=y(i)*fi(T,P,comp,K,y,zy,i)-x(i)*fi(T,P,comp,K,x,zx,i);
end
ansf(i+1)=1-sum(y);
answ=ansf;
% for i=1:size(x,2)
%     y(i)=x(i)*((fi(T,P,comp,K,x,zx,i)+eps)/(fi(T,P,comp,K,y,zy,i)+eps)+eps);
% end
% y
% answ=1-sum(y);
%answ=sum(y)-1;
% ansf(end+1)=sum(y)-1;

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

function H = Hv(comp,x,T)
for i=1:size(comp.Mw,2)
    Hi(i)=sum(cell2mat(comp.hv{i}).*[1 T^1 T^2 T^3 T^4 T^5])*comp.Mw(i)+comp.hoffset(i);%J/mol
end
H=sum(Hi.*x');

function H = Hl(comp,x,T,P)
R=8.314;
zv=1;%Zv(T,P,comp,K,y);
dT=1e-3;
for i=1:size(comp.Mw,2)
    a=cell2mat(comp.an{i}(1));b=cell2mat(comp.an{i}(2));c=cell2mat(comp.an{i}(3));
    d=cell2mat(comp.an{i}(4));e=cell2mat(comp.an{i}(5));f=cell2mat(comp.an{i}(6));
    P1=exp(a+b/(T+c)+d*log(T)+e*T^f)*1000;
    P2=exp(a+b/(T+dT+c)+d*log(T+dT)+e*(T+dT)^f)*1000;
    dpdT(i)=(P2-P1)/dT; %pa/k
end
dHlvi=dpdT*(R*T/P*zv)*T; %J/mol
dHlv=sum(dHlvi.*x');
H=(Hv(comp,x,T)-dHlv);