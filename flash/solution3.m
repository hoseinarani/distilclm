function x=solution3(J,e)
eror1=1;
eror2=0;
x=[0 0 0 0];
while abs(eror1-eror2)>1e-6
    eror2=eror1;
    x(1)=(-e(1)-(J(1,2)*x(2)+J(1,3)*x(3)+J(1,4)*x(4)))/J(1,1);
    x(2)=(-e(2)-(J(2,1)*x(1)+J(2,3)*x(3)+J(2,4)*x(4)))/J(2,2);
    x(3)=(-e(3)-(J(3,1)*x(1)+J(3,2)*x(2)+J(3,4)*x(4)))/J(3,3);
    x(4)=(-e(4)-(J(4,1)*x(1)+J(4,2)*x(2)+J(4,3)*x(3)))/J(4,4);
    eror1=norm(x);
end