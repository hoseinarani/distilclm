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