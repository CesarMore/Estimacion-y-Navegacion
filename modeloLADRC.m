
function[dxdt] = modeloLADRC(t,x,U1,U2,U3)

l=0.1;
kax=0.008;
kay=0.008;
Jz=0.1; 
kaz=0.009;
Jx=0.0552;
Jy=0.0552;
Jr=1;
ku=0.1188;
ks=0.0036;

    q1=100;
    q2=100;
    q3=100;
    q4=100;

d1 = 5;
d2 = 5;
d3 = 5;



dxdt(1,1) = x(2);
dxdt(2,1) = (-kay*x(2))/Jy + x(4)*Jr*(q1+q3-q2-q4)+(x(4)*x(6)*(Jz-Jx))/Jy+d2+(l*ku*U2)/Jy;
dxdt(3,1) = x(4);
dxdt(4,1) = (-kax*x(4))/Jx + x(2)*Jr*(q1+q3-q2-q4)+(x(2)*x(6)*(Jy-Jz))/Jx+d1+(l*ku*U1)/Jx;
dxdt(5,1) = x(6); 
dxdt(6,1) = (-kaz*x(6))/Jz +(x(4)*x(2)*(Jx-Jy))/Jz+d3+(ks*U3)/Jz;

end