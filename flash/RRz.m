function Rz=RRz(By,Bz,ky,kz,n)
Rz=sum(n.*(kz-1)./(1+By.*(ky-1)+Bz.*(kz-1)));