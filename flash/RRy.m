function Ry=RRy(By,Bz,ky,kz,n)
Ry=sum(n.*(ky-1)./(1+By.*(ky-1)+Bz.*(kz-1)));