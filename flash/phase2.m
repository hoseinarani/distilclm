function [Bx,By,Bz,x,y,z,zx,zy,zz]=phase2
clc
T = 290; %K
P = 101.3e3; %Pa
% C1, C2, C3, nC4, nC5, H2O, CO2, N2
index = [61,100,132,181,223,20,46,29];
n=[0.588235,0.058824,0.058824,0.058824,0.058824,0.058824,0.058824,0.058824];
%                       c2      c3  N-hexane water  xylene 
%[Bx,By,Bz,x,y,z,zx,zy,zz]=phase2(340,101.3e3,[    100    132    271    20   323     ],[.2 .2 .2 .2 .2])

comp.index=index;
comp.Mw=[datas(comp.index,1)]';
comp.Tfp=[datas(comp.index,2)]';
comp.Tb=[datas(comp.index,3)]';
comp.Tc=[datas(comp.index,4)]';
comp.Pc=[datas(comp.index,5)]'.*101.3e3;
comp.Zc=[datas(comp.index,7)]';
comp.Vc=[datas(comp.index,6)]';
comp.w=[datas(comp.index,8)]';
comp.CCpig=[datas(comp.index,12:15)]'.*4.184;
comp.miol=[datas(comp.index,16:17)]';
comp.dipol=[datas(comp.index,11)]';

K=zeros(size(n,2),size(n,2));
K=zeros(size(n,2),size(n,2));
R=8.314;
kz=zeros(size(n))+1e+6;
kz=comp.Pc./P.*exp(5.37.*(1+comp.w).*(1-comp.Tc./T));
[Bx,x,z,zx,zz]=twophaseVL(T,P,comp,K,n,kz);








