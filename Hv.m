x=ones(size(comp.Mw,2),1)/size(comp.Mw,2);
Hv(comp.hv,x,T,comp.Mw)
T=300
function H = Hv(comp,x,T)
for i=1:size(comp.Mw,2)
    Hi(i)=sum(cell2mat(comp.hv{i}).*[1 T^1 T^2 T^3 T^4 T^5])*x(i)*comp.Mw(i);%J/mol
end
H=sum(Hi);

function H = Hl(comp,x,T,P)
R=8.314;
zv=1;%Zv(T,P,comp,K,y);
dT=1e-5;
for i=1:size(comp.Mw,2)
    a=cell2mat(comp.an{i}(1));b=cell2mat(comp.an{i}(2));c=cell2mat(comp.an{i}(3));
    d=cell2mat(comp.an{i}(4));e=cell2mat(comp.an{i}(5));f=cell2mat(comp.an{i}(6));
    P1=exp(a+b/(T+c)+d*log(T)+e*T^f)*1000;
    P2=exp(a+b/(T+dT+c)+d*log(T+dT)+e*(T+dT)^f)*1000;
    dpdT(i)=(P2-P1)/dT; %pa/k
end
dHlv=sum(dpdT.*x)*(R*T/P*zv)*T; %J/mol
H=Hv(comp,x,T)-dHlv;